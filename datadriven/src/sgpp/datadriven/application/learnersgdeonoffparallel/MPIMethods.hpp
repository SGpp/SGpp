// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/PendingMPIRequest.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>
#include <list>

/**
 * Constant for determining whether the grid is currently receiving segmented update lists
 */
#define GRID_TEMPORARILY_INCONSISTENT 5
/**
 * Constant for determining whether receiving the list of deleted grid points has completed
 */
#define GRID_RECEIVED_DELETED_INDEXES 6
/**
 * Constant for determining whether receiving the list of added grid points has completed
 */
#define GRID_RECEIVED_ADDED_POINTS 7

/**
 * Check to determine whether a size_t to integer cast is safe
 */
#define CHECK_SIZE_T_TO_INT(x) if ((x) > INT_MAX)\
    {\
        throw sgpp::base::algorithm_exception("size_t to integer cast error");\
    }

/**
 * Check to determine whether a int to unsigned int cast is safe
 */
#define CHECK_INT_TO_UINT(x) if ((x) < 0)\
    {\
        throw sgpp::base::algorithm_exception("size_t to integer cast error");\
    }

namespace sgpp {
namespace datadriven {

// Forward declare Learner, as we use only pointer
class LearnerSGDEOnOffParallel;


/**
 * Structure that holds data for pending tracking requests
 * that need to be checked against incoming messages.
 */
struct MessageTrackRequest {
  /**
   * A predicate to apply to the incoming MPI request.
   * If predicate returns true, will increase the targetHits by one.
   */
  std::function<bool(PendingMPIRequest &)> predicate;
  /**
   * Counter to determine how many times the predicate has evaluated to true.
   */
  unsigned int targetHits;
  /**
   * Parameter to determine how many times the predicate should evaluate to true
   * before the track request completes.
   */
  unsigned int currentHits;
};

class MPIMethods {
 public:
  /**
   * Initializes the MPI setup for the specified learner.
   * This may only be called once, as it calls MPI_Init.
   * @param learnerInstance The instance of the learner to handle MPI communication for.
   */
  static void initMPI(LearnerSGDEOnOffParallel *learnerInstance);

  /**
   * Check whether the current role of this node is master
   * @return Whether the current role of this node is master
   */
  static bool isMaster();

  /**
   * Process any asynchronous MPI requests that completed but have not been processed yet.
   */
  static void processCompletedMPIRequests();

  /**
   * Process a specific completed MPI request.
   * @param pendingMPIRequest The completed MPI request.
   */
  static void processIncomingMPICommands(PendingMPIRequest &pendingMPIRequest);

  /**
   * Receive potentially segmented message that contains modifications to the grid.
   * Extract the data into a refinement result and apply to grid or apply to system matrix.
   * This includes additions, deletions, and system matrix updates.
   * @param networkMessage The message containing the modification data.
   */
  static void
  receiveGridComponentsUpdate(sgpp::datadriven::RefinementResultNetworkMessage *networkMessage);

  /**
   * Cancel any remaining requests and shutdown the MPI communicator.
   * May be called only once during a programs runtime.
   * Collectives (for example MPI_Ibcast) cannot be cancelled. These must complete before
   * calling this method.
   */
  static void finalizeMPI();

  /**
   * Broadcast a command that does not require any additional information
   * @param commandId The command to broadcast
   */
  static void bcastCommandNoArgs(MPI_COMMAND_ID commandId);

  /**
   * Tell a worker to learn from a batch.
   * @param workerID The MPI rank of the worker to command.
   * @param batchOffset The offset from the start of the training set to learn from.
   * @param batchSize The size of the batch to learn from.
   * @param doCrossValidation Whether to apply cross-validation
   */
  static void assignBatch(int workerID,
                          size_t batchOffset,
                          size_t batchSize,
                          bool doCrossValidation);

  /**
   * Get the number of nodes participating in the current problem
   * @return The number of participating nodes
   */
  static unsigned int getWorldSize();

  /**
   * Wait for any of the currently pending MPI requests to complete and process them.
   */
  static void waitForAnyMPIRequestsToComplete();

  /**
   * Send the results of training back to the master node for merging.
   * @param classIndex The index of the class to send results for.
   * @param batchOffset The offset from the start of the dataset used in this batch.
   * @param batchSize The size of the current batch.
   * @param alphaVector The results vector to transmit.
   * @return The number of successfully sent values.
   */
  static size_t sendMergeGridNetworkMessage(size_t classIndex, size_t batchOffset, size_t batchSize,
                                            base::DataVector &alphaVector);

  /**
   * Check how many pending MPI requests have been registered and not completed yet.
   * @return The number of pending MPI requests.
   */
  static size_t getQueueSize();

  /**
   * Copy the level index data for as many grid points as possible into the buffer specified.
   * @param buffer The start of the buffer to copy to.
   * @param bufferEnd The end of the buffer to copy to.
   * @param iterator The start of the data to copy.
   * @param listEnd The end of the data to copy.
   * @return The number of grid points copied.
   */
  static size_t fillBufferWithLevelIndexData(
      void *buffer, const void *bufferEnd,
      std::list<std::vector<LevelIndexPair>>::iterator &iterator,
      std::list<std::vector<LevelIndexPair>>::const_iterator &listEnd);

  /**
   * Send an assembled packet as an MPI Broadcast from the master.
   * @param mpiPacket The packet to transmit.
   * @return The now pending MPI transmission request.
   */
  static PendingMPIRequest &sendIBcast(MPI_Packet *mpiPacket);

  /**
   * Copy data from a specified iterator type into a message's payload.
   * @tparam Iterator The type of iterator used to read values to copy.
   * @param buffer The start of the buffer to copy to.
   * @param bufferEnd The end of the buffer to copy to.
   * @param iterator The start of the data to copy from.
   * @param listEnd The end of the data to copy from.
   * @return How many values of the iterator's type could be copied.
   */
  template<typename Iterator>
  static size_t fillBufferWithData(void *buffer, void *bufferEnd, Iterator &iterator,
                                   Iterator &listEnd);

  /**
   * Send the refinement updates for a specified class to all workers over a broadcast.
   * @param classIndex The index of the modified class.
   * @param deletedGridPointsIndices The list containing deleted grid points indices.
   * @param addedGridPoints The coordinates of any added grid points.
   */
  static void sendRefinementUpdates(size_t &classIndex, std::list<size_t> &deletedGridPointsIndices,
                                    std::list<LevelIndexVector> &addedGridPoints);

  /**
   * Send an MPI command using point to point communication without any further payload.
   * @param destinationRank The MPI rank of the receiving node.
   * @param commandId The command to send.
   */
  static void sendCommandNoArgs(int destinationRank, MPI_COMMAND_ID commandId);

  /**
   * Send a new system matrix decomposition with either broadcast or normal send.
   * @param classIndex The class index of the new system matrix decomposition
   * @param newSystemMatrixDecomposition The new system matrix decomposition
   * @param mpiTarget The MPI rank to send the decomposition to. Use MPI_ANY_SOURCE for broadcast.
   */
  static void
  sendSystemMatrixDecomposition(const size_t &classIndex,
                                sgpp::base::DataMatrix &newSystemMatrixDecomposition,
                                int mpiTarget);

  /**
   * Assign the task of updating the system matrix to the specified worker for processing.
   * @param workerID The MPI rank of the worker to assign to.
   * @param classIndex The class index of the system matrix to update.
   */
  static void assignSystemMatrixUpdate(int workerID, size_t classIndex);

  /**
   * Wait for an incoming message of the specified type.
   * Can optionally include a number of occurences to wait for
   * and a predicate to test each message against.
   * @param commandId The command id of the target packet.
   * @param numOccurrences The number of specified messages to receive.
   * @param predicate A predicate to test each message against.
   */
  static void waitForIncomingMessageType(MPI_COMMAND_ID commandId,
                                         unsigned int numOccurrences = 1,
                                         std::function<bool(
                                             PendingMPIRequest &)> predicate = [](
                                             PendingMPIRequest &request) { return true; });

  /**
   * Wait for a specific grid to return to its consistent state before continuing.
   * @param classIndex The grid's class index.
   */
  static void waitForGridConsistent(size_t classIndex);

  /**
   * Check whether there are any outgoing open MPI requests that have not completed yet
   * @return Whether such a request exists
   */
  static bool hasPendingOutgoingRequests();

  /**
    * Send an MPI packet to a target destination asynchronously.
    * @param destinationRank The MPI rank of the destination.
    * @param mpiPacket The MPI packet to send.
    * @param packetSize The size of the MPI packet to send, defaults to the entire packet.
    * @param highPriority Whether to send the packet on the high priority no_wait channel.
    * @return The pending MPI request to track completion.
    */
  static PendingMPIRequest &
  sendISend(int destinationRank, MPI_Packet *mpiPacket,
            size_t packetSize = sizeof(MPI_Packet),
            bool highPriority = false);

 protected:
  /**
   * Structure to track any pending tracking requests that need to be checked on every incoming
   * message.
   */
  static std::list<MessageTrackRequest> messageTrackRequests;

  /**
   * Structure to track all pending MPI requests that have been registered but not completed.
   */
  static std::list<PendingMPIRequest> pendingMPIRequests;

  /**
   * Structure to hold the actual MPI_Request structures in memory sequentially
   * along with the additional information for each request.
   */
  static MPIRequestPool mpiRequestStorage;

  /**
   * The number of participating nodes for solving the problem.
   */
  static unsigned int mpiWorldSize;

  /**
   * The instance of the learner for which MPI communications are currently being handled.
   */
  static LearnerSGDEOnOffParallel *learnerInstance;

  /**
   * Learn from a batch based on the instructions in the specified packet.
   * @param assignBatchMessage The packet specifying learning parameters.
   */
  static void runBatch(MPI_Packet *assignBatchMessage);

  /**
   * Receive a packet of changes from the master and apply them to the grid and refinement result
   * or to the system matrix decomposition.
   * @param networkMessage The packet containing the changes
   */
  static void
  receiveMergeGridNetworkMessage(MergeGridNetworkMessage &networkMessage);

  /**
   * Create a pending MPI request for the specified packet in preparation to sending it.
   * @param mpiPacket The packet of data for the request.
   * @param isInbound Whether the request is an inbound request.
   * @return Reference to the newly created PendingMPIRequest.
   */
  static PendingMPIRequest &createPendingMPIRequest(MPI_Packet *mpiPacket, bool isInbound);

  /**
   * Process a specific completed MPI request from its position in the request list.
   * @param pendingMPIRequestIterator The position of the completed MPI request.
   */
  static void
  processCompletedMPIRequest(
      const std::list<sgpp::datadriven::PendingMPIRequest>::iterator &pendingMPIRequestIterator);

  /**
   * Find a PendingMPIRequest from the index of the completed MPI_Request structure in
   * the pool storage.
   * @param completedRequestIndex The index of the completed MPI_Request in the MPIRequestPool.
   * @return The located PendingMPIRequest.
   */
  static std::list<sgpp::datadriven::PendingMPIRequest>::iterator findPendingMPIRequest(
      unsigned int completedRequestIndex);

  /**
   * Wait for any of the requests in the MPI request pool to complete.
   */
  static unsigned int executeMPIWaitAny();

  /**
   * Callback function that processes the request and zeros then memory region afterwards
   */
  static void handleIncomingRequestFromCallback(PendingMPIRequest &request);

  /**
   * Create a message track request that tests each message against a predicate until the target
   * number of occurences are found.
   * @param numOccurrences How many successful matches are required.
   * @param predicate The predicate to test messages against.
   * @return The newly allocated track request's position in the storage list
   */
  static std::list<sgpp::datadriven::MessageTrackRequest>::iterator
  createTrackRequest(unsigned int numOccurrences,
                     const std::function<bool(PendingMPIRequest &)> &predicate);

  /**
   * Calculate the size of an MPI_Packet based on the size of its contained payload.
   * @param containedPacketSize The size of the contained payload.
   * @return The total MPI_Packet size.
   */
  static size_t calculateTotalPacketSize(size_t containedPacketSize);
};
}  // namespace datadriven
}  // namespace sgpp
