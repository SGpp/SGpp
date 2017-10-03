//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_MPIMETHODS_H
#define SGPP_MPIMETHODS_H

namespace sgpp {
    namespace datadriven {
        class MPIMethods {
        public:
            static void initMPI();

            static void synchronizeEndOfDataPass();

            static bool isMaster();

            static void sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults);


            static void waitForMPIRequestsToComplete();

            static void processCompletedMPIRequests();

            static void processIncomingMPICommands(LearnerSGDEOnOffParallel *learnerInstance, MPI_Packet *mpiPacket);


            static void
            receiveGridComponentsUpdate(LearnerSGDEOnOffParallel *learnerInstance,
                                        RefinementResultNetworkMessage *networkMessage);

        protected:
            //Pending MPI Requests
            static std::vector<PendingMPIRequest> pendingMPIRequests;
            static int mpiWorldSize;

            static void startSynchronizingPackets();

            void initMPI(LearnerSGDEOnOffParallel *learnerInstance);

            template<typename DataType>
            static size_t fillBufferWithData(void *buffer, void *bufferEnd, std::list::iterator iterator,
                                             std::list::iterator listEnd);

            template<typename DataType>
            static size_t
            fillBufferWithVectorData(void *buffer, void *bufferEnd, std::list<std::vector>::iterator iterator,
                                     std::list::iterator listEnd);

            static void sendRefinementUpdates(size_t classIndex, const RefinementResultsUpdateType &updateType,
                                              const std::list::iterator &iterator,
                                              const std::list<unsigned long>::iterator &listEnd);

            static size_t sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
                                                     const RefinementResult &refinementResult, int offset,
                                                     std::list::iterator &iterator);

            void sendMergeGridCommand(std::vector<base::DataVector *> &alphas);
        };
    }
}


#endif //SGPP_MPIMETHODS_H
