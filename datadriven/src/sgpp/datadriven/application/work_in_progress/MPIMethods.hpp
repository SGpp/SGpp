//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_MPIMETHODS_H
#define SGPP_MPIMETHODS_H

#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>

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

            void initMPI(LearnerSGDEOnOffParallel *learnerInstance);

            template<typename DataType>
            static size_t fillBufferWithData(void *buffer, void *bufferEnd, typename DataType::const_iterator iterator,
                                             typename DataType::const_iterator listEnd);

            template<typename DataType>
            static size_t
            fillBufferWithVectorData(void *buffer, void *bufferEnd,
                                     typename std::vector<DataType>::const_iterator iterator,
                                     typename DataType::const_iterator listEnd);

            template<typename DataType>
            static void sendRefinementUpdates(size_t &classIndex, const RefinementResultsUpdateType updateType,
                                              typename std::list<DataType>::const_iterator &iterator,
                                              typename std::list<DataType>::const_iterator &listEnd);

//            static size_t sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
//                                                     const RefinementResult &refinementResult, int offset,
//                                                     std::list::iterator &iterator);

            void sendMergeGridCommand(std::vector<base::DataVector *> &alphas);
        };
    }
}


#endif //SGPP_MPIMETHODS_H
