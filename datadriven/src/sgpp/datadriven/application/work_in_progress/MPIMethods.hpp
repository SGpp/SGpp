//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_MPIMETHODS_H
#define SGPP_MPIMETHODS_H

#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>
#include "PendingMPIRequest.hpp"

#define GRID_TEMPORARILY_INCONSISTENT 5
#define GRID_RECEIVED_DELETED_INDEXES 6
#define GRID_RECEIVED_ADDED_POINTS 7

#define CHECK_SIZE_T_TO_INT(x) if((x) > INT_MAX)\
    {\
        throw sgpp::base::algorithm_exception("size_t to integer cast error");\
    }

#define CHECK_INT_TO_UINT(x) if((x) < 0)\
    {\
        throw sgpp::base::algorithm_exception("size_t to integer cast error");\
    }

namespace sgpp {
    namespace datadriven {
        struct MessageTrackRequest {
            std::function<bool(PendingMPIRequest &)> predicate;
            unsigned int targetHits;
            unsigned int currentHits;
        };

        class MPIMethods {
        public:

            static void initMPI(LearnerSGDEOnOffParallel *learnerInstance);

            static bool isMaster();

//            static void sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults);

            static void processCompletedMPIRequests();

            static void processIncomingMPICommands(PendingMPIRequest &pendingMPIRequest);


            static void
            receiveGridComponentsUpdate(sgpp::datadriven::RefinementResultNetworkMessage *networkMessage);

            static void finalizeMPI();

            static void bcastCommandNoArgs(MPI_COMMAND_ID commandId);

            static void assignBatch(int workerID, size_t batchOffset, size_t batchSize, bool doCrossValidation);

        static unsigned int getWorldSize();

            static void waitForAnyMPIRequestsToComplete();

            static size_t sendMergeGridNetworkMessage(size_t classIndex, size_t batchOffset, size_t batchSize,
                                                      base::DataVector &alphaVector);

            static size_t getQueueSize();

            static size_t fillBufferWithLevelIndexData(void *buffer, const void *bufferEnd,
                                                       std::list<std::vector<LevelIndexPair>>::iterator &iterator,
                                                       std::list<std::vector<LevelIndexPair>>::const_iterator &listEnd);

            static PendingMPIRequest & sendIBcast(MPI_Packet *mpiPacket);

            template<typename Iterator>
            static size_t fillBufferWithData(void *buffer, void *bufferEnd, Iterator &iterator,
                                             Iterator &listEnd);

            static void sendRefinementUpdates(size_t &classIndex, std::list<size_t> &deletedGridPointsIndexes,
                                              std::list<LevelIndexVector> &addedGridPoints);


            static void sendCommandNoArgs(int destinationRank, MPI_COMMAND_ID commandId);

            static void
            sendCholeskyDecomposition(const size_t &classIndex, DataMatrix &newCholeskyDecomposition, int mpiTarget);

        static void assignSystemMatrixUpdate(int workerID, size_t classIndex);

        static void waitForIncomingMessageType(MPI_COMMAND_ID commandId,
                                               unsigned int numOccurrences = 1,
                                               std::function<bool(
                                                       PendingMPIRequest &)> predicate = [](
                                                           PendingMPIRequest &request) { return true; });

            static void waitForGridConsistent(size_t classIndex);

        protected:

            static std::list<MessageTrackRequest> messageTrackRequests;

            //Pending MPI Requests
            static std::list<PendingMPIRequest> pendingMPIRequests;
            static MPIRequestPool mpiRequestStorage;
        static unsigned int mpiWorldSize;
            static LearnerSGDEOnOffParallel *learnerInstance;


//            static size_t sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
//                                                     const RefinementResult &refinementResult, int offset,
//                                                     std::list::iterator &iterator);

            static PendingMPIRequest &
            sendISend(int destinationRank, MPI_Packet *mpiPacket,
                      size_t packetSize = sizeof(MPI_Packet),
                      bool highPriority = false);

            static void runBatch(MPI_Packet *pPacket);

            static void
            receiveMergeGridNetworkMessage(MergeGridNetworkMessage &networkMessage);

            static PendingMPIRequest &createPendingMPIRequest(MPI_Packet *mpiPacket, bool isInbound);

            static void
            processCompletedMPIRequest(
                    const std::list<sgpp::datadriven::PendingMPIRequest>::iterator &pendingMPIRequestIterator);

            static std::list<sgpp::datadriven::PendingMPIRequest>::iterator findPendingMPIRequest(
                    unsigned int completedRequestIndex);

        static unsigned int executeMPIWaitAny();

            static void handleIncomingRequestFromCallback(PendingMPIRequest &request);

            static std::list<sgpp::datadriven::MessageTrackRequest>::iterator
            createTrackRequest(unsigned int numOccurrences, const std::function<bool(PendingMPIRequest &)> &predicate);

            static size_t calculateTotalPacketSize(size_t containedPacketSize);
        };
    }
}


#endif //SGPP_MPIMETHODS_H
