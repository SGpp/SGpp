// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <mpi.h>

#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <thread>
#include <climits>
#include <cstring>
#include <list>
#include <numeric>

namespace sgpp {
namespace datadriven {

// Pending MPI Requests
std::list<PendingMPIRequest> MPIMethods::pendingMPIRequests;
MPIRequestPool MPIMethods::mpiRequestStorage;
unsigned int MPIMethods::mpiWorldSize = 0;
LearnerSGDEOnOffParallel *MPIMethods::learnerInstance;
std::list<MessageTrackRequest> MPIMethods::messageTrackRequests;

bool MPIMethods::isMaster() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == MPI_MASTER_RANK;
}

void MPIMethods::initMPI(LearnerSGDEOnOffParallel *learnerInstance) {
  MPI_Init(nullptr, nullptr);

  int signedMPIWorldSize = -1;

  // Get World Size
  MPI_Comm_size(MPI_COMM_WORLD, &signedMPIWorldSize);

  CHECK_INT_TO_UINT(signedMPIWorldSize);
  mpiWorldSize = static_cast<unsigned int>(signedMPIWorldSize);

  // Get Rank
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get Processor Name
  char mpiProcessorName[MPI_MAX_PROCESSOR_NAME_LENGTH];
  int nameLength;
  MPI_Get_processor_name(mpiProcessorName, &nameLength);

  std::cout << "Processor " << mpiProcessorName << " (rank " << world_rank
            << ") has joined MPI pool of size " << mpiWorldSize << std::endl;

  // Setup receiving messages from master/workers
  {
    auto *mpiPacket = new MPI_Packet;
    auto &unicastInputRequest = createPendingMPIRequest(mpiPacket, true);
    unicastInputRequest.disposeAfterCallback = false;
    unicastInputRequest.callback = [](PendingMPIRequest &request) {
      D(std::cout << "Incoming MPI unicast" << std::endl;)
      handleIncomingRequestFromCallback(request);

      D(std::cout << "Restarting irecv request." << std::endl;)
      MPI_Irecv(request.buffer, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
                MPI_TAG_STANDARD_COMMAND, MPI_COMM_WORLD, request.getMPIRequestFromHandle());
    };

    MPI_Irecv(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
              MPI_TAG_STANDARD_COMMAND,
              MPI_COMM_WORLD,
              unicastInputRequest.getMPIRequestFromHandle());

    std::cout << "Started listening for unicasts from any sources" << std::endl;
  }
  // Setup receiving high priority messages from master/workers
  {
    auto *mpiPacket = new MPI_Packet;
    auto &unicastInputRequest = createPendingMPIRequest(mpiPacket, true);
    unicastInputRequest.disposeAfterCallback = false;
    unicastInputRequest.callback = [](PendingMPIRequest &request) {
      D(std::cout << "Incoming MPI high priority unicast" << std::endl;)
      handleIncomingRequestFromCallback(request);

      D(std::cout << "Restarting irecv request." << std::endl;)
      MPI_Irecv(request.buffer, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
                MPI_TAG_HIGH_PRIORITY_NO_BLOCK, MPI_COMM_WORLD,
                request.getMPIRequestFromHandle());
    };

    MPI_Irecv(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
              MPI_TAG_HIGH_PRIORITY_NO_BLOCK,
              MPI_COMM_WORLD,
              unicastInputRequest.getMPIRequestFromHandle());

    std::cout << "Started listening for high priority unicasts from any sources" << std::endl;
  }
  if (!isMaster()) {
    auto *mpiPacket = new MPI_Packet;
    auto &broadcastInputRequest = createPendingMPIRequest(mpiPacket, true);
    broadcastInputRequest.disposeAfterCallback = false;
    broadcastInputRequest.callback = [](PendingMPIRequest &request) {
      D(std::cout << "Incoming MPI broadcast" << std::endl;)
      handleIncomingRequestFromCallback(request);

      D(std::cout << "Restarting ibcast request." << std::endl;)
      MPI_Ibcast(request.buffer, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                 MPI_COMM_WORLD, request.getMPIRequestFromHandle());
    };
    MPI_Ibcast(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
               MPI_COMM_WORLD,
               broadcastInputRequest.getMPIRequestFromHandle());

    std::cout << "Started listening for broadcasts from task master" << std::endl;
  }

  MPIMethods::learnerInstance = learnerInstance;
}

void MPIMethods::handleIncomingRequestFromCallback(PendingMPIRequest &request) {
  processIncomingMPICommands(request);

  D(std::cout << "Zeroing MPI Request" << std::endl;)
//            std::memset(request.getMPIRequestFromHandle(), 0, sizeof(MPI_Request));

  D(std::cout << "Zeroing Buffer" << std::endl;)
  std::memset(request.buffer, 0, sizeof(MPI_Packet));
}

void
MPIMethods::sendRefinementUpdates(size_t &classIndex, std::list<size_t> &deletedGridPointsIndices,
                                  std::list<LevelIndexVector> &addedGridPoints) {
  // Deleted grid points
  {
    std::list<size_t>::const_iterator iterator = deletedGridPointsIndices.begin();
    std::list<size_t>::const_iterator listEnd = deletedGridPointsIndices.end();
    while (iterator != listEnd) {
      auto *mpiPacket = new MPI_Packet;
      mpiPacket->commandID = UPDATE_GRID;

      auto *networkMessage =
          static_cast<RefinementResultNetworkMessage *>(static_cast<void *>(mpiPacket->payload));

      networkMessage->classIndex = classIndex;
      networkMessage->updateType = DELETED_GRID_POINTS_LIST;

      size_t numPointsInBuffer = fillBufferWithData<std::list<size_t>::const_iterator>(
          static_cast<void *>(networkMessage->payload),
          static_cast<void *>(std::end(networkMessage->payload)),
          iterator,
          listEnd);
      networkMessage->listLength = numPointsInBuffer;

      networkMessage->gridversion = (iterator == listEnd) ? GRID_RECEIVED_DELETED_INDEXES
                                                          : GRID_TEMPORARILY_INCONSISTENT;

      D(std::cout << "Sending updated for class " << networkMessage->classIndex
                  << " with " << networkMessage->listLength
                  << " deletions" << " (grid version " << networkMessage->gridversion << ")"
                  << std::endl;)

      sendIBcast(mpiPacket);
    }
  }
  // Added grid points
  {
    auto iterator = addedGridPoints.begin();
    std::list<LevelIndexVector>::const_iterator listEnd = addedGridPoints.end();
    while (iterator != listEnd) {
      auto *mpiPacket = new MPI_Packet;
      mpiPacket->commandID = UPDATE_GRID;

      auto *networkMessage =
          static_cast<RefinementResultNetworkMessage *>(static_cast<void *>(mpiPacket->payload));

      networkMessage->classIndex = classIndex;
      networkMessage->updateType = ADDED_GRID_POINTS_LIST;

      size_t numPointsInBuffer = fillBufferWithLevelIndexData(networkMessage->payload,
                                                              std::end(
                                                                  networkMessage->payload),
                                                              iterator,
                                                              listEnd);
      networkMessage->listLength = numPointsInBuffer;

      networkMessage->gridversion = (iterator == listEnd) ? GRID_RECEIVED_ADDED_POINTS
                                                          : GRID_TEMPORARILY_INCONSISTENT;


      D(std::cout << "Sending updated for class " << networkMessage->classIndex
                  << " with " << networkMessage->listLength
                  << " additions" << " (grid version " << networkMessage->gridversion << ")"
                  << std::endl;)

      sendIBcast(mpiPacket);
    }
  }
}

// USE MPI_ANY_SOURCE to send a broadcast from master
void MPIMethods::sendSystemMatrixDecomposition(const size_t &classIndex,
                                               DataMatrix &newSystemMatrixDecomposition,
                                               int mpiTarget) {
  auto iterator = newSystemMatrixDecomposition.begin();
  auto listEnd = newSystemMatrixDecomposition.end();

  size_t offset = 0;
  while (iterator != listEnd) {
    auto *mpiPacket = new MPI_Packet;
    mpiPacket->commandID = UPDATE_GRID;

    auto *networkMessage =
        static_cast<RefinementResultNetworkMessage *>(static_cast<void *>(mpiPacket->payload));

    networkMessage->classIndex = classIndex;
    networkMessage->updateType = SYSTEM_MATRIX_DECOMPOSITION;

    auto *systemMatrixNetworkMessage =
        static_cast<RefinementResultSystemMatrixNetworkMessage *>(
            static_cast<void *>(networkMessage->payload));

    systemMatrixNetworkMessage->matrixWidth = newSystemMatrixDecomposition.getNcols();
    systemMatrixNetworkMessage->matrixHeight = newSystemMatrixDecomposition.getNrows();
    systemMatrixNetworkMessage->offset = offset;

    size_t numPointInBuffer = fillBufferWithData(systemMatrixNetworkMessage->payload,
                                                 std::end(systemMatrixNetworkMessage->payload),
                                                 iterator, listEnd);
    networkMessage->listLength = numPointInBuffer;

    offset += numPointInBuffer;

    networkMessage->gridversion = (iterator == listEnd) ? learnerInstance->getLocalGridVersion(
        classIndex) : GRID_TEMPORARILY_INCONSISTENT;

    if (mpiTarget != MPI_ANY_SOURCE) {
      D(std::cout << "Sending system matrix for class " << networkMessage->classIndex
                  << " offset " << systemMatrixNetworkMessage->offset
                  << " with " << networkMessage->listLength
                  << " values" << " (grid version " << networkMessage->gridversion
                  << ", target "
                  << mpiTarget << ")" << std::endl;)
      sendISend(mpiTarget, mpiPacket, sizeof(MPI_Packet), true);
    } else {
      D(std::cout << "Broadcasting system matrix for class " << networkMessage->classIndex
                  << " offset " << systemMatrixNetworkMessage->offset
                  << " with " << networkMessage->listLength
                  << " values" << " (grid version " << networkMessage->gridversion << ")"
                  << std::endl;)
      sendIBcast(mpiPacket);
    }
  }
}

void MPIMethods::bcastCommandNoArgs(MPI_COMMAND_ID commandId) {
  auto *mpiPacket = new MPI_Packet;
  mpiPacket->commandID = commandId;

  sendIBcast(mpiPacket);
}

void MPIMethods::sendCommandNoArgs(const int destinationRank, MPI_COMMAND_ID commandId) {
  auto *mpiPacket = new MPI_Packet;
  mpiPacket->commandID = commandId;

  sendISend(destinationRank, mpiPacket);
}

PendingMPIRequest &MPIMethods::sendIBcast(MPI_Packet *mpiPacket) {
  PendingMPIRequest &pendingMPIRequest = createPendingMPIRequest(mpiPacket, false);

  MPI_Ibcast(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
             MPI_COMM_WORLD, pendingMPIRequest.getMPIRequestFromHandle());

  D(std::cout << "Ibcast request stored at " << &pendingMPIRequest << std::endl;)
  return pendingMPIRequest;
}

PendingMPIRequest &
MPIMethods::sendISend(const int destinationRank, MPI_Packet *mpiPacket, size_t packetSize,
                      bool highPriority) {
  PendingMPIRequest &pendingMPIRequest = createPendingMPIRequest(mpiPacket, false);

  CHECK_SIZE_T_TO_INT(packetSize)
  // Point to the request in vector instead of stack
  MPI_Isend(mpiPacket, static_cast<int >(packetSize), MPI_UNSIGNED_CHAR, destinationRank,
            highPriority ? MPI_TAG_HIGH_PRIORITY_NO_BLOCK : MPI_TAG_STANDARD_COMMAND,
            MPI_COMM_WORLD, pendingMPIRequest.getMPIRequestFromHandle());

  D(std::cout << "Isend request stored at " << &pendingMPIRequest << std::endl;)
  return pendingMPIRequest;
}

PendingMPIRequest &MPIMethods::createPendingMPIRequest(MPI_Packet *mpiPacket, bool isInbound) {
  pendingMPIRequests.emplace_back(&mpiRequestStorage);
  PendingMPIRequest &pendingMPIRequest = pendingMPIRequests.back();
  pendingMPIRequest.disposeAfterCallback = true;
  pendingMPIRequest.inbound = isInbound;
  pendingMPIRequest.callback = [](PendingMPIRequest &request) {
    D(std::cout << "Pending MPI request " << &request << " ("
                << (request.inbound ? "inbound" : "outbound")
                << ") completed." << std::endl;)
  };
  pendingMPIRequest.buffer = mpiPacket;
  return pendingMPIRequest;
}

void MPIMethods::receiveMergeGridNetworkMessage(MergeGridNetworkMessage &networkMessage) {
  base::DataVector alphaVector(networkMessage.alphaTotalSize);

  auto *payload = static_cast<double *>(static_cast<void *>(networkMessage.payload));
  for (size_t index = 0; index < networkMessage.payloadLength; index++) {
    alphaVector[networkMessage.payloadOffset + index] = payload[index];
  }

  bool isLastPacketInSeries =
      networkMessage.payloadLength + networkMessage.payloadOffset ==
          networkMessage.alphaTotalSize;

  learnerInstance->mergeAlphaValues(networkMessage.classIndex, networkMessage.gridversion,
                                    alphaVector,
                                    networkMessage.batchOffset, networkMessage.batchSize,
                                    isLastPacketInSeries);

  D(std::cout << "Updated alpha values from network message offset "
              << networkMessage.payloadOffset
              << ", class " << networkMessage.classIndex
              << ", length " << networkMessage.payloadLength << ", alpha vector length "
              << networkMessage.alphaTotalSize << std::endl;)
}

size_t
MPIMethods::sendMergeGridNetworkMessage(size_t classIndex, size_t batchOffset, size_t batchSize,
                                        base::DataVector &alphaVector) {
  size_t offset = 0;
  auto beginIterator = alphaVector.begin();
  auto endIterator = alphaVector.end();

  while (offset < alphaVector.size()) {
    auto *mpiPacket = new MPI_Packet;
    mpiPacket->commandID = MERGE_GRID;

    void *payloadPointer = &(mpiPacket->payload);

    auto *networkMessage = static_cast<MergeGridNetworkMessage *>(payloadPointer);

    networkMessage->classIndex = classIndex;
    networkMessage->gridversion = learnerInstance->getLocalGridVersion(classIndex);
    networkMessage->payloadOffset = offset;
    networkMessage->batchSize = batchSize;
    networkMessage->batchOffset = batchOffset;
    networkMessage->alphaTotalSize = alphaVector.size();

    size_t numPointsInPacket = 0;

    numPointsInPacket = fillBufferWithData(networkMessage->payload,
                                           std::end(networkMessage->payload),
                                           beginIterator, endIterator);

    networkMessage->payloadLength = numPointsInPacket;

    D(
        std::cout << "Sending merge for class " << classIndex
                  << " offset " << offset
                  << " with " << numPointsInPacket << " values"
                  << " and grid version " << networkMessage->gridversion << std::endl;
        std::cout << "Alpha sum is "
                  << std::accumulate(alphaVector.begin(), alphaVector.end(), 0.0)
                  << std::endl;
    )
    sendISend(MPI_MASTER_RANK, mpiPacket);
    offset += numPointsInPacket;
  }

  return offset;
}

size_t
MPIMethods::fillBufferWithLevelIndexData(void *buffer, const void *bufferEnd,
                                         std::list<LevelIndexVector>::iterator &iterator,
                                         std::list<LevelIndexVector>::const_iterator &listEnd) {
  size_t copiedValues = 0;
  auto *bufferPointer = static_cast<LevelIndexPair *>(buffer);
  while (iterator != listEnd && bufferPointer < bufferEnd) {
    LevelIndexVector levelIndexVector = *iterator;
    if (bufferPointer + 2 * levelIndexVector.size() * sizeof(size_t) >= bufferEnd) {
      break;
    }
    for (auto &levelIndexPair : levelIndexVector) {
      bufferPointer->level = levelIndexPair.level;
      bufferPointer->index = levelIndexPair.index;
      bufferPointer++;
    }
    copiedValues++;
    iterator++;
  }
  return copiedValues;
}

template<typename Iterator>
size_t
MPIMethods::fillBufferWithData(void *buffer, void *bufferEnd, Iterator &iterator,
                               Iterator &listEnd) {
  auto *bufferPointer = (typename std::iterator_traits<Iterator>::value_type *) buffer;
  size_t copiedValues = 0;
  while (bufferPointer + sizeof(typename std::iterator_traits<Iterator>::value_type) <
      bufferEnd &&
      iterator != listEnd) {
    *bufferPointer = *iterator;

    bufferPointer++;
    iterator++;
    copiedValues++;
  }
  return copiedValues;
}

void MPIMethods::processCompletedMPIRequests() {
  D(std::cout << "Checking " << pendingMPIRequests.size() << " pending MPI requests"
              << std::endl;)
  auto pendingMPIRequestIterator = pendingMPIRequests.end();
  auto listBegin = pendingMPIRequests.begin();

  // In order to process the send requests first we start from the back

  while (pendingMPIRequestIterator != listBegin) {
    pendingMPIRequestIterator--;

    MPI_Status mpiStatus{};
    int operationCompleted;

    D(std::cout << "Testing request " << &*pendingMPIRequestIterator << std::endl;)
    if (MPI_Test(pendingMPIRequestIterator->getMPIRequestFromHandle(), &operationCompleted,
                 &mpiStatus) !=
        MPI_SUCCESS) {
      std::cout << "Error MPI Test reported" << std::endl;
      throw sgpp::base::algorithm_exception("MPI Test returned error.");
    }
    D(std::cout << "Pending request has status " << operationCompleted
                << " MPI_ERROR: " << mpiStatus.MPI_ERROR
                << " MPI SOURCE: " << mpiStatus.MPI_SOURCE
                << " MPI TAG: " << mpiStatus.MPI_TAG << std::endl;)
    if (operationCompleted != 0) {
      processCompletedMPIRequest(pendingMPIRequestIterator);


      D(std::cout << "Relaunching processCompletedMPIRequests" << std::endl;)
      processCompletedMPIRequests();
      break;
    }
  }
}

void MPIMethods::processCompletedMPIRequest(
    const std::list<sgpp::datadriven::PendingMPIRequest>::iterator &pendingMPIRequestIterator) {
  D(std::cout << "Updating " << messageTrackRequests.size() << " track requests" << std::endl;)
  for (MessageTrackRequest &trackRequest : messageTrackRequests) {
    if (trackRequest.predicate(*pendingMPIRequestIterator)) {
      trackRequest.currentHits++;
    }
  }

  D(std::cout << "Executing callback" << std::endl;)
  // Execute the callback
  pendingMPIRequestIterator->callback(*pendingMPIRequestIterator);
  D(std::cout << "Callback complete" << std::endl;)


  if (pendingMPIRequestIterator->disposeAfterCallback) {
    D(std::cout << "Attempting to delete pending mpi request" << std::endl;)
    delete pendingMPIRequestIterator->buffer;

    pendingMPIRequests.erase(pendingMPIRequestIterator);
    D(std::cout << "Deleted pending mpi request" << std::endl;)
  }
}

void MPIMethods::waitForAnyMPIRequestsToComplete() {
  unsigned int completedRequest = executeMPIWaitAny();
  processCompletedMPIRequest(findPendingMPIRequest(completedRequest));
}

unsigned int MPIMethods::executeMPIWaitAny() {
  int completedRequest = -1;
  MPI_Status mpiStatus{};
  D(std::cout << "Waiting for " << pendingMPIRequests.size() << " MPI requests to complete"
              << std::endl;)

  D(auto begin = std::chrono::high_resolution_clock::now();)

  CHECK_SIZE_T_TO_INT(mpiRequestStorage.size());
  int requestStatus = MPI_Waitany(static_cast<int>(mpiRequestStorage.size()),
                                  mpiRequestStorage.getMPIRequests(),
                                  &completedRequest,
                                  &mpiStatus);
  D(
      auto end = std::chrono::high_resolution_clock::now();

      auto deltaMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
          end - begin).count();
      if (deltaMilliseconds > 100) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            begin.time_since_epoch()).count() - 1
                  << " active 1" << std::endl;
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            begin.time_since_epoch()).count()
                  << " active 0" << std::endl;
        std::cout << "Idle time " << deltaMilliseconds << "ms" << std::endl;
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            end.time_since_epoch()).count()
                  << " active 0" << std::endl;
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            end.time_since_epoch()).count() + 1
                  << " active 1" << std::endl;
      }
  )
  D(std::cout << "MPI request " << completedRequest << " completed"
              << " MPI_ERROR: " << mpiStatus.MPI_ERROR
              << " MPI SOURCE: " << mpiStatus.MPI_SOURCE
              << " MPI TAG: " << mpiStatus.MPI_TAG << " MPI REQUEST STATUS: " << requestStatus
              << std::endl;)

  CHECK_SIZE_T_TO_INT(mpiRequestStorage.size());
  if (requestStatus != MPI_SUCCESS
      || completedRequest < 0 ||
      completedRequest > static_cast<int>(mpiRequestStorage.size())) {
    std::cout << "Error: Completed requests returned invalid index " << completedRequest
              << " operation status " << requestStatus << std::endl;
    throw sgpp::base::algorithm_exception("MPI Request completed without a valid index");
  }
  CHECK_INT_TO_UINT(completedRequest)
  return static_cast<unsigned int>(completedRequest);
}

void MPIMethods::waitForIncomingMessageType(MPI_COMMAND_ID commandId, unsigned int numOccurrences,
                                            std::function<bool(PendingMPIRequest &)> predicate) {
  D(std::cout << "Waiting for " << numOccurrences << " messages of type " << commandId
              << std::endl;)
  auto trackRequest = createTrackRequest(numOccurrences, [commandId, predicate](
      PendingMPIRequest &mpiRequest) {
    return mpiRequest.inbound && mpiRequest.buffer->commandID == commandId &&
        predicate(mpiRequest);
  });

  while (trackRequest->currentHits < trackRequest->targetHits) {
    unsigned int completedRequest = executeMPIWaitAny();
    const auto &pendingMPIRequestIterator = findPendingMPIRequest(completedRequest);

    processCompletedMPIRequest(pendingMPIRequestIterator);
    D(std::cout << "Received " << trackRequest->currentHits << "/" << trackRequest->targetHits
                << " messages of type " << commandId
                << std::endl;)
  }

  D(std::cout << "Received all " << trackRequest->targetHits << " messages of type " << commandId
              << std::endl;)
  messageTrackRequests.erase(trackRequest);
}

std::list<sgpp::datadriven::MessageTrackRequest>::iterator
MPIMethods::createTrackRequest(unsigned int numOccurrences,
                               const std::function<bool(PendingMPIRequest &)> &predicate) {
  messageTrackRequests.emplace_front();
  auto trackRequest = messageTrackRequests.begin();
  trackRequest->currentHits = 0;
  trackRequest->predicate = predicate;
  trackRequest->targetHits = numOccurrences;
  return trackRequest;
}

std::list<PendingMPIRequest>::iterator
MPIMethods::findPendingMPIRequest(unsigned int completedRequestIndex) {
  auto iterator = pendingMPIRequests.begin();
  const std::list<sgpp::datadriven::PendingMPIRequest>::iterator
      &listEnd = pendingMPIRequests.end();
  while (iterator != listEnd) {
    if (iterator->getMPIRequestIndex() == completedRequestIndex) {
      return iterator;
    }
    iterator++;
  }
  std::cout << "Pending MPI Request not found: " << completedRequestIndex << std::endl;
  throw sgpp::base::algorithm_exception("Could not find PendingMPIRequest for Pool Index.");
}

void MPIMethods::receiveGridComponentsUpdate(RefinementResultNetworkMessage *networkMessage) {
  size_t classIndex = networkMessage->classIndex;
  RefinementResult &refinementResult = learnerInstance->getRefinementHandler().getRefinementResult(
      classIndex);

  size_t listLength = networkMessage->listLength;
  size_t processedPoints = 0;
  void *bufferEnd = std::end(networkMessage->payload);

  D(std::cout << "Receiving " << listLength << " grid modifications for class " << classIndex
              << " (update type " << networkMessage->updateType << ", remote grid version "
              << networkMessage->gridversion << ")"
              << std::endl;)

  if (learnerInstance->checkGridStateConsistent(classIndex) &&
      !learnerInstance->isVersionConsistent(networkMessage->gridversion)) {
    D(std::cout << "Received first message in multi segment grid update of type "
                << networkMessage->updateType << std::endl;)
    while (!isMaster() && networkMessage->updateType != SYSTEM_MATRIX_DECOMPOSITION
        && (!refinementResult.addedGridPoints.empty() ||
            !refinementResult.deletedGridPointsIndices.empty())) {
      std::cout
          << "Received first message in multi segment grid update,"
          << " however refinement results have not been cleared ("
          << refinementResult.addedGridPoints.size() << " additions, "
          << refinementResult.deletedGridPointsIndices.size() << " deletions)."
          << std::endl;
      waitForIncomingMessageType(ASSIGN_BATCH);
    }
  }
  switch (networkMessage->updateType) {
    case DELETED_GRID_POINTS_LIST: {
      auto *bufferIterator = static_cast<size_t *>(static_cast<void *>(networkMessage->payload));
      while (bufferIterator < bufferEnd && processedPoints < listLength) {
        refinementResult.deletedGridPointsIndices.push_back(*bufferIterator);
        bufferIterator++;
        processedPoints++;
      }
    }
      break;
    case ADDED_GRID_POINTS_LIST: {
      auto *bufferIterator =
          static_cast<LevelIndexPair *>(static_cast<void *>(networkMessage->payload));
      size_t dimensionality = learnerInstance->getDimensionality();
      while (bufferIterator < bufferEnd && processedPoints < listLength) {
        LevelIndexVector dataVector(dimensionality);
        size_t currentDimension = 0;
        while (currentDimension < dimensionality && bufferIterator < bufferEnd) {
          dataVector[currentDimension].level = bufferIterator->level;
          dataVector[currentDimension].index = bufferIterator->index;
          bufferIterator++;

          currentDimension++;
        }
        refinementResult.addedGridPoints.push_back(dataVector);
        processedPoints++;
      }
      break;
    }
    case SYSTEM_MATRIX_DECOMPOSITION: {
      auto *systemMatrixNetworkMessage =
          static_cast<RefinementResultSystemMatrixNetworkMessage *>(
              static_cast<void *>(networkMessage->payload));
      DataMatrix &systemMatrixDecomposition =
          learnerInstance->getDensityFunctions()[classIndex]
              .first->getOfflineObject().getDecomposedMatrix();

      D(std::cout << "Receiving system matrix decomposition update at offset "
                  << systemMatrixNetworkMessage->offset
                  << std::endl;)

      if (systemMatrixNetworkMessage->offset == 0) {
        size_t oldSize = systemMatrixDecomposition.size();
        systemMatrixDecomposition.resizeRowsCols(systemMatrixNetworkMessage->matrixHeight,
                                             systemMatrixNetworkMessage->matrixWidth);
        std::cout << "Adjusted size of system matrix decomposition " << classIndex << " from "
                  << oldSize
                  << " to "
                  << systemMatrixDecomposition.size() << std::endl;
      } else if (learnerInstance->getLocalGridVersion(classIndex) !=
          GRID_TEMPORARILY_INCONSISTENT) {
        std::cout
            << "Update with non-null offset arrived "
            << "but grid version not set to inconsistent (version is "
            << learnerInstance->getLocalGridVersion(classIndex) << std::endl;
        throw sgpp::base::algorithm_exception("Secondary update received on consistent grid.");
      }

      auto *bufferIterator =
          static_cast<double *>(static_cast<void *>(systemMatrixNetworkMessage->payload));
      while (bufferIterator < bufferEnd && processedPoints < listLength) {
//                        std::cout << "Copy from " << &*bufferIterator << " to "
//                                  << (systemMatrixNetworkMessage->offset + processedPoints)
// << std::endl;
        systemMatrixDecomposition[systemMatrixNetworkMessage->offset +
            processedPoints] = *bufferIterator;
        bufferIterator++;
        processedPoints++;
      }

//                    std::cout << "Copy successful" << std::endl;
      if (isMaster() &&
          systemMatrixNetworkMessage->offset + processedPoints ==
              systemMatrixDecomposition.size()) {
        learnerInstance->setLocalGridVersion(classIndex, networkMessage->gridversion);

        D(std::cout << "Received system matrix decomposition for class " << classIndex
                    << ", will now broadcast decomposition" << std::endl;)
        sendSystemMatrixDecomposition(classIndex, systemMatrixDecomposition, MPI_ANY_SOURCE);
      }
      break;
    }
    default: {
      std::cout << "Received an update request with unknown id " << networkMessage->updateType
                << std::endl;
      throw sgpp::base::algorithm_exception("Update request with unknown ID received.");
    }
  }
  D(std::cout << "Updated refinement result or system matrix decomposition " << classIndex << " ("
              << refinementResult.addedGridPoints.size()
              << " additions, "
              << refinementResult.deletedGridPointsIndices.size() <<
              " deletions)" << std::endl;)

  // If this is not the last message in a series (gridversion inconsistent),
  // then don't update variables yet
  if (networkMessage->gridversion == GRID_RECEIVED_ADDED_POINTS && !isMaster()) {
    D(std::cout << "Updating class " << classIndex
                << " variables as grid is now consistent with version "
                << networkMessage->gridversion << std::endl;)
    learnerInstance->getRefinementHandler()
        .updateClassVariablesAfterRefinement(classIndex,
                                             &refinementResult,
                                             learnerInstance->getDensityFunctions()[classIndex]
                                                 .first.get(),
                                                 learnerInstance->getGrid(classIndex));
  }
  learnerInstance->setLocalGridVersion(classIndex, networkMessage->gridversion);
}

void MPIMethods::processIncomingMPICommands(PendingMPIRequest &pendingMPIRequest) {
  MPI_Packet *mpiPacket = pendingMPIRequest.buffer;
  D(std::cout << "Processing incoming command " << mpiPacket->commandID << std::endl;)
  void *networkMessagePointer = &(mpiPacket->payload);

  switch (mpiPacket->commandID) {
    case UPDATE_GRID: {
      auto *refinementResultNetworkMessage =
          static_cast<RefinementResultNetworkMessage *>(networkMessagePointer);
      receiveGridComponentsUpdate(refinementResultNetworkMessage);
    }
      break;
    case MERGE_GRID: {
      auto *mergeGridNetworkMessage = static_cast<MergeGridNetworkMessage *>(networkMessagePointer);
      receiveMergeGridNetworkMessage(*mergeGridNetworkMessage);
    }
      break;
    case ASSIGN_BATCH:runBatch(mpiPacket);
      break;
    case COMPUTE_UPDATE_SYSTEM_MATRIX_DECOMPOSITION: {
      auto *message = static_cast<AssignSystemMatrixUpdateNetworkMessage *>(networkMessagePointer);
      learnerInstance->computeNewSystemMatrixDecomposition(message->classIndex,
                                                       message->gridversion);
      break;
    }
    case SHUTDOWN:std::cout << "Worker shutdown requested" << std::endl;
      learnerInstance->shutdownMPINodes();
      std::cout << "Marking pending mpi request for dispose" << std::endl;
      pendingMPIRequest.disposeAfterCallback = true;
      break;
    case WORKER_SHUTDOWN_SUCCESS:std::cout << "Worker has acknowledged shutdown" << std::endl;
      break;
    case nullptr_COMMAND:std::cout << "Error: Incoming command has undefined command id" << std::endl;
      throw sgpp::base::algorithm_exception("MPI_Packet with nullptr command received.");
    default:std::cout << "Error: MPI unknown command id: " << mpiPacket->commandID << std::endl;
      throw sgpp::base::algorithm_exception("MPI_Packet with unknown command received.");
  }
}

void MPIMethods::finalizeMPI() {
  std::cout << pendingMPIRequests.size() << " MPI requests pending before finalize" << std::endl;
  size_t requestNum = 0;
  for (PendingMPIRequest &pendingMPIRequest : pendingMPIRequests) {
    std::cout << "Cancelling pending mpi request " << requestNum << " at " << &pendingMPIRequest
              << std::endl;
    MPI_Request *mpiRequestHandle = pendingMPIRequest.getMPIRequestFromHandle();
    MPI_Cancel(mpiRequestHandle);
    MPI_Wait(mpiRequestHandle, MPI_STATUS_IGNORE);
    delete pendingMPIRequest.buffer;
    requestNum++;
  }
  pendingMPIRequests.clear();
  std::cout << "Finalizing MPI" << std::endl;
  MPI_Finalize();
}

void MPIMethods::assignBatch(const int workerID, size_t batchOffset, size_t batchSize,
                             bool doCrossValidation) {
  auto *mpiPacket = new MPI_Packet;
  mpiPacket->commandID = ASSIGN_BATCH;

  auto *message = static_cast<AssignBatchNetworkMessage *>(static_cast<void *>(mpiPacket->payload));
  message->batchOffset = batchOffset;
  message->batchSize = batchSize;
  message->doCrossValidation = doCrossValidation;

  sendISend(workerID, mpiPacket, calculateTotalPacketSize(sizeof(AssignBatchNetworkMessage)));
}

size_t MPIMethods::calculateTotalPacketSize(size_t containedPacketSize) {
  return offsetof(MPI_Packet, MPI_Packet::payload) + containedPacketSize;
}

void MPIMethods::runBatch(MPI_Packet *assignBatchMessage) {
  auto *message =
      static_cast<AssignBatchNetworkMessage *>(static_cast<void *>(assignBatchMessage->payload));
  D(std::cout << "runbatch dim " << learnerInstance->getDimensionality() << std::endl;)
  D(std::cout << "creating dataset" << std::endl;)
  Dataset dataset{message->batchSize, learnerInstance->getDimensionality()};
  learnerInstance->workBatch(dataset, message->batchOffset, message->doCrossValidation);
}

void MPIMethods::assignSystemMatrixUpdate(const int workerID, size_t classIndex) {
  auto *mpiPacket = new MPI_Packet;
  mpiPacket->commandID = COMPUTE_UPDATE_SYSTEM_MATRIX_DECOMPOSITION;

  auto *message =
      static_cast<AssignSystemMatrixUpdateNetworkMessage *>(
          static_cast<void *>(mpiPacket->payload));
  message->gridversion = learnerInstance->getLocalGridVersion(classIndex);
  message->classIndex = classIndex;

  sendISend(workerID, mpiPacket,
            calculateTotalPacketSize(sizeof(AssignSystemMatrixUpdateNetworkMessage)),
            true);
}

unsigned int MPIMethods::getWorldSize() {
  return mpiWorldSize;
}

size_t MPIMethods::getQueueSize() {
  return pendingMPIRequests.size();
}

void MPIMethods::waitForGridConsistent(size_t classIndex) {
  while (!learnerInstance->checkGridStateConsistent(classIndex)) {
    D(std::cout << "Grid " << classIndex << " is not yet consistent (version "
                << learnerInstance->getLocalGridVersion(classIndex) << ") waiting for update."
                << std::endl;)
    waitForIncomingMessageType(UPDATE_GRID);
  }
}

bool MPIMethods::hasPendingOutgoingRequests() {
  return std::any_of(pendingMPIRequests.begin(),
                     pendingMPIRequests.end(),
                     [](PendingMPIRequest &request) { return !request.inbound; });
}

}  // namespace datadriven
}  // namespace sgpp
