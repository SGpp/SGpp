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

            static void synchronizeEndOfDataPass();

            static bool isMaster();

            static void sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults);


            static void waitForMPIRequestsToComplete();

            static void processCompletedMPIRequests();

            static void processIncomingMPICommands(LearnerSGDEOnOffParallel *learnerInstance,
                                                   sgpp::datadriven::MPI_Packet *mpiPacket);


            static void
            receiveGridComponentsUpdate(LearnerSGDEOnOffParallel *learnerInstance,
                                        sgpp::datadriven::RefinementResultNetworkMessage *networkMessage);

        protected:
            //Pending MPI Requests
            static std::vector<sgpp::datadriven::PendingMPIRequest> pendingMPIRequests;
            static int mpiWorldSize;

            static void startSynchronizingPackets();


            template<typename Iterator>
            static size_t fillBufferWithData(void *buffer, void *bufferEnd, Iterator iterator,
                                             Iterator listEnd);

            template<typename Iterator, typename ValueType>
            static size_t
            fillBufferWithVectorData(void *buffer, const void *bufferEnd,
                                     Iterator iterator,
                                     Iterator listEnd, size_t sizeOfDataType);

            template<typename Iterator>
            static void sendRefinementUpdates(size_t &classIndex, const RefinementResultsUpdateType updateType,
                                              Iterator &iterator,
                                              Iterator &listEnd);

//            static size_t sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
//                                                     const RefinementResult &refinementResult, int offset,
//                                                     std::list::iterator &iterator);

        };
    }
}


#endif //SGPP_MPIMETHODS_H
