//
// Created by Vincent_Bode on 11.09.2017.
//

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp>

#define BOOST_TEST_DYN_LINK
#define SCHEDULER_BATCH_SIZE 1337
#define TEST_DIMENSION 1
#define TEST_DATASET_SIZE 1

BOOST_AUTO_TEST_SUITE(MPIMethods_Test)

using namespace sgpp::datadriven;

sgpp::datadriven::LearnerSGDEOnOffParallel *learnerInstance;

void freeInstance() {
  free(learnerInstance);
}

void createInstance() {
  if (learnerInstance == nullptr) {
    sgpp::datadriven::DBMatDensityConfiguration dbMatDensityConf{};
    sgpp::datadriven::Dataset trainData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::datadriven::Dataset testData(TEST_DATASET_SIZE, TEST_DIMENSION);
    sgpp::base::DataVector classLabels(2);
    classLabels[0] = -1;
    classLabels[1] = 1;
    RoundRobinScheduler roundRobinScheduler = RoundRobinScheduler(SCHEDULER_BATCH_SIZE);
    learnerInstance = new sgpp::datadriven::LearnerSGDEOnOffParallel(dbMatDensityConf,
                                                                     trainData, testData, nullptr,
                                                                     classLabels, 2, false,
                                                                     0.0, 0.01,
                                                                     roundRobinScheduler
    );
    atexit(freeInstance);
  }
}

BOOST_AUTO_TEST_CASE(AssignBatchTest) {
  createInstance();

  BOOST_REQUIRE(MPIMethods::isMaster());
  BOOST_REQUIRE(MPIMethods::getWorldSize() == 1);

  size_t batchOffset = 330;
  bool doCrossValidation = false;
  learnerInstance->assignBatchToWorker(batchOffset, doCrossValidation);

  MPIMethods::waitForIncomingMessageType(ASSIGN_BATCH,
                                         1,
                                         [batchOffset, doCrossValidation](PendingMPIRequest &pendingMPIRequest) {
                                           auto message =
                                               static_cast<AssignBatchNetworkMessage *>(static_cast<void *>(pendingMPIRequest.buffer->payload));
                                           return message->batchOffset == batchOffset
                                               && message->doCrossValidation == doCrossValidation
                                               && message->batchSize == SCHEDULER_BATCH_SIZE;
                                         });
}

BOOST_AUTO_TEST_CASE(MergeAlphaValuesTest) {
  createInstance();

  BOOST_REQUIRE(MPIMethods::isMaster());
  BOOST_REQUIRE(MPIMethods::getWorldSize() == 1);

  MPI_Packet mpiPacket{};
  mpiPacket.commandID = MERGE_GRID;

  auto *message = static_cast<MergeGridNetworkMessage *>(static_cast<void *>(mpiPacket.payload));

  message->batchSize = SCHEDULER_BATCH_SIZE;
  message->batchOffset = 50;
  message->gridversion = 10;
  message->classIndex = 0;
  message->batchSize = 1;
  message->alphaTotalSize = 2;
  message->payloadOffset = 1;
  message->payloadLength = 1;
}

BOOST_AUTO_TEST_SUITE_END()
