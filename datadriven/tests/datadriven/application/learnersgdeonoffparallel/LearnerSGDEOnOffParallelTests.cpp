// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI
#ifndef USE_SCALAPACK  // this test interferes with the ScaLAPACK tests

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>

#define SCHEDULER_BATCH_SIZE 1337
#define TEST_DIMENSION 3
#define TEST_DATASET_SIZE 3

BOOST_AUTO_TEST_SUITE(MPIMethods_Test)

using sgpp::datadriven::AssignTaskResult;
using sgpp::datadriven::DataMatrix;
using sgpp::datadriven::LearnerSGDEOnOffParallel;
using sgpp::datadriven::LevelIndexPair;
using sgpp::datadriven::LevelIndexVector;
using sgpp::datadriven::MERGE_GRID;
using sgpp::datadriven::MergeGridNetworkMessage;
using sgpp::datadriven::MPI_Packet;
using sgpp::datadriven::MPIMethods;
using sgpp::datadriven::MPIRequestPool;
using sgpp::datadriven::RefinementResult;
using sgpp::datadriven::RefinementResultNetworkMessage;
using sgpp::datadriven::RefinementResultSystemMatrixNetworkMessage;
using sgpp::datadriven::RoundRobinScheduler;
using sgpp::datadriven::TRAIN_FROM_BATCH;
using sgpp::datadriven::UPDATE_GRID;

sgpp::datadriven::LearnerSGDEOnOffParallel *learnerInstance;
sgpp::datadriven::RoundRobinScheduler *scheduler;

/**
 * Ensure that specific messages are not too long to be wrapped in the parent packet.
 * These are used to prevent errors by adding variables to messages without shortening
 * their payload.
 */
static_assert(sizeof(MergeGridNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Merge Grid Network Message too long.");
static_assert(sizeof(RefinementResultNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Refinement result Network Message too long.");
static_assert(sizeof(RefinementResultSystemMatrixNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Refinement result Cholesky Network Message too long.");

void freeInstance() {
  delete learnerInstance;
  delete scheduler;
}

void createInstance() {
  if (learnerInstance == nullptr) {
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = TEST_DIMENSION;
    gridConfig.level_ = 3;
    gridConfig.type_ = sgpp::base::GridType::Linear;

    sgpp::datadriven::RegularizationConfiguration regularizationConfig{};
    regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
    regularizationConfig.lambda_ = 0.01;

    sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
    densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::DenseIchol;

    sgpp::base::AdaptivityConfiguration adaptivityConfig;
    adaptivityConfig.numRefinements_ = 2;
    adaptivityConfig.numRefinementPoints_ = 7;
    adaptivityConfig.refinementThreshold_ = 0.0;  // only required for surplus refinement

    sgpp::datadriven::Dataset trainData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::datadriven::Dataset testData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::base::DataVector classLabels(2);
    classLabels[0] = -1;
    classLabels[1] = 1;
    scheduler = new RoundRobinScheduler(SCHEDULER_BATCH_SIZE);
    learnerInstance = new sgpp::datadriven::LearnerSGDEOnOffParallel(
        gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig, trainData,
        testData, nullptr, classLabels, 2, false, 0.0, *scheduler);

    atexit(freeInstance);
  }
}

void sendMergeGridPacket(size_t batchSize, size_t batchOffset, size_t gridversion,
                         size_t classIndex, size_t alphaTotalSize, size_t payloadOffset,
                         size_t payloadLength) {
  auto *mpiPacket = new MPI_Packet;
  mpiPacket->commandID = MERGE_GRID;

  auto *message = static_cast<MergeGridNetworkMessage *>(static_cast<void *>(mpiPacket->payload));

  message->batchOffset = batchOffset;
  message->gridversion = gridversion;
  message->classIndex = classIndex;
  message->batchSize = batchSize;
  message->alphaTotalSize = alphaTotalSize;
  message->payloadOffset = payloadOffset;
  message->payloadLength = payloadLength;

  MPIMethods::sendISend(0, mpiPacket);
}

BOOST_AUTO_TEST_CASE(RoundRobinSchedulerTest) {
  createInstance();

  BOOST_CHECK(scheduler->isReadyForRefinement());

  AssignTaskResult result{};
  scheduler->assignTaskVariableTaskSize(TRAIN_FROM_BATCH, result);

  BOOST_CHECK(result.workerID == 1);
  BOOST_CHECK(result.taskSize == SCHEDULER_BATCH_SIZE);

  BOOST_CHECK(scheduler->isReadyForRefinement());
  scheduler->onRefinementStarted();

  BOOST_CHECK(!scheduler->isReadyForRefinement());
  BOOST_CHECK_THROW(scheduler->onRefinementStarted(), sgpp::base::algorithm_exception);

  BOOST_CHECK_THROW(scheduler->onMergeRequestIncoming(0, SCHEDULER_BATCH_SIZE, 10, 12),
                    sgpp::base::algorithm_exception);

  // Once for each class
  scheduler->onMergeRequestIncoming(0, SCHEDULER_BATCH_SIZE, 10, 11);
  BOOST_CHECK(!scheduler->isReadyForRefinement());
  scheduler->onMergeRequestIncoming(0, SCHEDULER_BATCH_SIZE, 10, 11);
  BOOST_CHECK(scheduler->isReadyForRefinement());
}

BOOST_AUTO_TEST_CASE(AssignBatchTest) {
  createInstance();

  BOOST_REQUIRE(MPIMethods::isMaster());
  BOOST_REQUIRE(MPIMethods::getWorldSize() == 1);

  //  size_t batchOffset = 330;
  //  bool doCrossValidation = false;
  //  learnerInstance->assignBatchToWorker(batchOffset, doCrossValidation);
  //
  //  MPIMethods::waitForIncomingMessageType(ASSIGN_BATCH,
  //     1,
  //     [batchOffset, doCrossValidation](PendingMPIRequest &pendingMPIRequest) {
  //       auto message =
  //           static_cast<AssignBatchNetworkMessage *>
  //                (static_cast<void *>(pendingMPIRequest.buffer->payload));
  //       return message->batchOffset == batchOffset
  //           && message->doCrossValidation == doCrossValidation
  //           && message->batchSize == SCHEDULER_BATCH_SIZE;
  //     });
}

BOOST_AUTO_TEST_CASE(SendRefinementResultsTest) {
  createInstance();

  size_t classIndex = 0;

  RefinementResult &refinementResult =
      learnerInstance->getRefinementHandler().getRefinementResult(classIndex);

  // Test deleted grid points
  refinementResult.deletedGridPointsIndices.emplace_back(1);
  refinementResult.deletedGridPointsIndices.emplace_back(3);
  refinementResult.deletedGridPointsIndices.emplace_back(5);

  MPIMethods::sendRefinementUpdates(classIndex, refinementResult.deletedGridPointsIndices,
                                    refinementResult.addedGridPoints);
  BOOST_CHECK(MPIMethods::hasPendingOutgoingRequests());
  MPIMethods::waitForAnyMPIRequestsToComplete();

  refinementResult.deletedGridPointsIndices.clear();

  // Test added grid points
  LevelIndexVector levelIndexVector{};
  LevelIndexPair levelIndexPair{3, 3};
  levelIndexVector.emplace_back(levelIndexPair);
  refinementResult.addedGridPoints.emplace_back(levelIndexVector);

  MPIMethods::sendRefinementUpdates(classIndex, refinementResult.deletedGridPointsIndices,
                                    refinementResult.addedGridPoints);
  BOOST_CHECK(MPIMethods::hasPendingOutgoingRequests());
  MPIMethods::waitForAnyMPIRequestsToComplete();
}

BOOST_AUTO_TEST_CASE(SendSystemMatrixDecompositionTest) {
  size_t classIndex = 0;
  DataMatrix systemMatrix(2, 2, -1.0);

  MPIMethods::sendSystemMatrixDecomposition(classIndex, systemMatrix, 0);
  BOOST_CHECK(MPIMethods::hasPendingOutgoingRequests());
  MPIMethods::waitForIncomingMessageType(UPDATE_GRID);

  DataMatrix &installedMatrix = learnerInstance->getDensityFunctions()[classIndex]
                                    .first->getOfflineObject()
                                    .getDecomposedMatrix();
  BOOST_CHECK(installedMatrix.size() == systemMatrix.size());
  BOOST_CHECK(installedMatrix[0] == systemMatrix[0]);
}

BOOST_AUTO_TEST_CASE(MergeAlphaValuesTest) {
  createInstance();

  BOOST_REQUIRE(MPIMethods::isMaster());
  BOOST_REQUIRE(MPIMethods::getWorldSize() == 1);

  // Test for correct message current version
  learnerInstance->setLocalGridVersion(0, 10);
  sendMergeGridPacket(1, 50, 10, 0, 31, 0, 3);

  MPIMethods::waitForIncomingMessageType(MERGE_GRID);

  // Test for correct message last version, no ref data
  learnerInstance->setLocalGridVersion(0, 11);
  sendMergeGridPacket(1, 50, 10, 0, 31, 0, 3);

  BOOST_CHECK_THROW(MPIMethods::waitForIncomingMessageType(MERGE_GRID),
                    sgpp::base::algorithm_exception);
  // Clean up
  learnerInstance->setLocalGridVersion(0, 10);
  MPIMethods::processCompletedMPIRequests();

  // Test for correct message last version, with ref data
  learnerInstance->setLocalGridVersion(0, 11);
  LevelIndexVector levelIndexVector{};
  LevelIndexPair levelIndexPair{3, 3};
  levelIndexVector.emplace_back(levelIndexPair);
  RefinementResult &refinementResult =
      learnerInstance->getRefinementHandler().getRefinementResult(0);
  refinementResult.addedGridPoints.clear();
  refinementResult.deletedGridPointsIndices.clear();
  refinementResult.addedGridPoints.emplace_back(levelIndexVector);
  sendMergeGridPacket(1, 50, 10, 0, 30, 0, 3);

  MPIMethods::waitForIncomingMessageType(MERGE_GRID);
  refinementResult.addedGridPoints.clear();

  // Test for correct message version -2 should end in error
  learnerInstance->setLocalGridVersion(0, 12);
  sendMergeGridPacket(1, 50, 10, 0, 31, 0, 3);

  BOOST_CHECK_THROW(MPIMethods::waitForIncomingMessageType(MERGE_GRID),
                    sgpp::base::algorithm_exception);
  // Clean up
  learnerInstance->setLocalGridVersion(0, 10);
  MPIMethods::processCompletedMPIRequests();
}

BOOST_AUTO_TEST_CASE(MPIRequestPool_Test) {
  sgpp::datadriven::MPIRequestPool requestPool;

  BOOST_CHECK(requestPool.size() == 0);
  requestPool.createMPIRequestHandle();
  BOOST_CHECK(requestPool.size() == 1);
  requestPool.createMPIRequestHandle();
  BOOST_CHECK(requestPool.size() == 2);

  requestPool.deleteMPIRequestHandle(1);
  BOOST_CHECK(requestPool.size() == 1);

  requestPool.createMPIRequestHandle();
  requestPool.createMPIRequestHandle();
  BOOST_CHECK(requestPool.size() == 3);

  // This should create blank space in the storage but not modify the pool size
  requestPool.deleteMPIRequestHandle(0);
  BOOST_CHECK(requestPool.size() == 3);

  // This should refill blank space in the storage but not modify the pool size
  requestPool.createMPIRequestHandle();
  BOOST_CHECK(requestPool.size() == 3);

  // This should create blank space in the storage but not modify the pool size
  requestPool.deleteMPIRequestHandle(1);
  BOOST_CHECK(requestPool.size() == 3);

  // This should cause the entire pool to shrink
  requestPool.deleteMPIRequestHandle(2);
  BOOST_CHECK(requestPool.size() == 1);

  // The current implementation cannot shrink back to size 0
  requestPool.deleteMPIRequestHandle(0);
  BOOST_CHECK(requestPool.size() <= 1);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SCALAPACK */
#endif /* USE_MPI */
