// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#define BOOST_TEST_DYN_LINK
#define SCHEDULER_BATCH_SIZE 1337
#define TEST_DIMENSION 3
#define TEST_DATASET_SIZE 3

BOOST_AUTO_TEST_SUITE(MPIMethods_Test)

using namespace sgpp::datadriven;

sgpp::datadriven::LearnerSGDEOnOffParallel *learnerInstance;
sgpp::datadriven::RoundRobinScheduler *scheduler;

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
    regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

    sgpp::datadriven::DBMatDecompostionType
        dt = sgpp::datadriven::DBMatDecompostionType::DenseIchol;

    sgpp::base::AdpativityConfiguration adaptConfig;
    adaptConfig.numRefinements_ = 2;
    adaptConfig.noPoints_ = 7;
    adaptConfig.threshold_ = 0.0;  // only required for surplus refinement

    sgpp::datadriven::DBMatDensityConfiguration dbMatDensityConf(gridConfig,
                                                                 adaptConfig,
                                                                 regularizationConfig.regType_,
                                                                 0.01, dt);
    sgpp::datadriven::Dataset trainData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::datadriven::Dataset testData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::base::DataVector classLabels(2);
    classLabels[0] = -1;
    classLabels[1] = 1;
    scheduler = new RoundRobinScheduler(SCHEDULER_BATCH_SIZE);
    learnerInstance = new sgpp::datadriven::LearnerSGDEOnOffParallel(dbMatDensityConf,
                                                                     trainData, testData, nullptr,
                                                                     classLabels, 2, false,
                                                                     0.0, 0.01,
                                                                     *scheduler);

    size_t nextCvStep = 50000;
    double cvLambdaStart = 1e-1;
    double cvLambdaEnd = 1e-10;
    int cvLambdaSteps = 10;
    bool cvLogScale = true;
    sgpp::base::DataMatrix *cvTestData = &testData.getData();
    sgpp::base::DataMatrix *cvTestDataRes = nullptr;  // needed?
    learnerInstance->setCrossValidationParameters(cvLambdaSteps,
                                                  cvLambdaStart,
                                                  cvLambdaEnd,
                                                  cvTestData,
                                                  cvTestDataRes,
                                                  cvLogScale);

    atexit(freeInstance);
  }
}

void sendMergeGridPacket(size_t batchSize,
                         size_t batchOffset,
                         size_t gridversion,
                         size_t classIndex,
                         size_t alphaTotalSize,
                         size_t payloadOffset,
                         size_t payloadLength) {
  MPI_Packet *mpiPacket = new MPI_Packet;
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

  RefinementResult &refinementResult = learnerInstance->getRefinementHandler()
      .getRefinementResult(0);
  // TODO(bodevt): Resume here
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
  learnerInstance->getRefinementHandler().
      getRefinementResult(0).addedGridPoints.emplace_back(levelIndexVector);
  sendMergeGridPacket(1, 50, 10, 0, 30, 0, 3);

  MPIMethods::waitForIncomingMessageType(MERGE_GRID);
  learnerInstance->getRefinementHandler().
      getRefinementResult(0).addedGridPoints.clear();

  // Test for correct message version -2 should end in error
  learnerInstance->setLocalGridVersion(0, 12);
  sendMergeGridPacket(1, 50, 10, 0, 31, 0, 3);

  BOOST_CHECK_THROW(MPIMethods::waitForIncomingMessageType(MERGE_GRID),
                    sgpp::base::algorithm_exception);
  // Clean up
  learnerInstance->setLocalGridVersion(0, 10);
  MPIMethods::processCompletedMPIRequests();
}

BOOST_AUTO_TEST_SUITE_END()



