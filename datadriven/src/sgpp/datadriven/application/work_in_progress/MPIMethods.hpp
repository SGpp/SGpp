//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_MPIMETHODS_H
#define SGPP_MPIMETHODS_H

#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>

namespace sgpp {
    namespace datadriven {
        class MPIMethods {
        public:

            static void initMPI(LearnerSGDEOnOffParallel *learnerInstance);

            static bool isMaster();

            static void sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults);

            static void processCompletedMPIRequests();

            static void processIncomingMPICommands(sgpp::datadriven::MPI_Packet *mpiPacket);


            static void
            receiveGridComponentsUpdate(sgpp::datadriven::RefinementResultNetworkMessage *networkMessage);

            static void synchronizeBarrier();

            static void finalizeMPI();

            static void bcastCommandNoArgs(MPI_COMMAND_ID commandId);

            static void assignBatch(const int workerID, size_t batchOffset, size_t batchSize, bool doCrossValidation);

            static int getWorldSize();

            static void waitForAnyMPIRequestsToComplete();

            static size_t sendMergeGridNetworkMessage(size_t classIndex, size_t batchSize, DataVector &alphaVector);

        protected:
            //Pending MPI Requests
            static std::list<sgpp::datadriven::PendingMPIRequest> pendingMPIRequests;
            static int mpiWorldSize;
            static LearnerSGDEOnOffParallel *learnerInstance;

            static void startSynchronizingPackets();

            static void sendIBcast(MPI_Packet *mpiPacket);

            template<typename Iterator>
            static size_t fillBufferWithData(void *buffer, void *bufferEnd, Iterator &iterator,
                                             Iterator &listEnd);

            template<typename Iterator, typename ValueType>
            static size_t
            fillBufferWithVectorData(void *buffer, const void *bufferEnd,
                                     Iterator &iterator,
                                     Iterator &listEnd, size_t sizeOfDataType);

            template<typename Iterator>
            static void sendRefinementUpdates(size_t &classIndex, RefinementResultsUpdateType updateType,
                                              Iterator &iterator,
                                              Iterator &listEnd);

//            static size_t sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
//                                                     const RefinementResult &refinementResult, int offset,
//                                                     std::list::iterator &iterator);

            static void endSynchronizingPackets();

            static void sendISend(const int destinationRank, MPI_Packet *mpiPacket);

            void sendCommandNoArgs(const int destinationRank, MPI_COMMAND_ID commandId);

            static void runBatch(MPI_Packet *pPacket);

            static size_t
            receiveMergeGridNetworkMessage(MergeGridNetworkMessage &networkMessage);

            void waitForAllMPIRequestsToComplete();

            static PendingMPIRequest &createPendingMPIRequest(MPI_Packet *mpiPacket);

        };
    }
}


#endif //SGPP_MPIMETHODS_H
