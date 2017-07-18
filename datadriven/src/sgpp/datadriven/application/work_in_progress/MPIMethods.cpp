//
// Created by Vincent_Bode on 29.06.2017.
//

#include <ctime>
#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <cstring>
#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
    namespace datadriven {

        bool MPIMethods::isMaster() {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            return rank == MPI_MASTER_RANK;
        }

        void MPIMethods::initMPI(LearnerSGDEOnOffParallel *learnerInstance) {
            MPI_Init(NULL, NULL);

            // Get World Size
            MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

            //Get Rank
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            //Get Processor Name
            char mpiProcessorName[MPI_MAX_PROCESSOR_NAME_LENGTH];
            int nameLength;
            MPI_Get_processor_name(mpiProcessorName, &nameLength);

            printf("Processor %s (rank %i) has joined MPI pool of size %i\n", mpiProcessorName, world_rank,
                   mpiWorldSize);

            //Setup receiving messages from master/workers
            {
                sgpp::datadriven::PendingMPIRequest unicastInputRequest;
                MPI_Packet *mpiPacket = new MPI_Packet;
                unicastInputRequest.buffer = mpiPacket;
                unicastInputRequest.disposeAfterCallback = false;
                unicastInputRequest.callback = [&learnerInstance](PendingMPIRequest &request) {
                    processIncomingMPICommands(learnerInstance, request.buffer);
                    MPI_Irecv(request.buffer, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
                              MPI_ANY_TAG, MPI_COMM_WORLD, &(request.request));
                };

                MPI_Irecv(mpiPacket, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                          MPI_COMM_WORLD,
                          &(unicastInputRequest.request));

                pendingMPIRequests.push_back(unicastInputRequest);
                printf("Started listening for unicasts from any sources\n");
            }
            if (!isMaster()) {
                PendingMPIRequest broadcastInputRequest;
                MPI_Packet *mpiPacket = new MPI_Packet;
                broadcastInputRequest.buffer = mpiPacket;
                broadcastInputRequest.disposeAfterCallback = false;
                broadcastInputRequest.callback = [&learnerInstance](PendingMPIRequest &request) {
                    processIncomingMPICommands(learnerInstance, request.buffer);
                    MPI_Ibcast(request.buffer, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                               MPI_COMM_WORLD, &(request.request));
                };
                MPI_Ibcast(mpiPacket, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK, MPI_COMM_WORLD,
                           &(broadcastInputRequest.request));
                pendingMPIRequests.push_back(broadcastInputRequest);

                printf("Started listening for broadcasts from task master\n");
            }

        }

        void MPIMethods::synchronizeEndOfDataPass() {
            throw sgpp::base::application_exception("Not implemented");
//            sgpp::parallel::myGlobalMPIComm->Barrier();
        }

        void MPIMethods::sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults) {

            for (size_t classIndex = 0; classIndex < refinementResults->size(); classIndex++) {
                RefinementResult refinementResult = (*refinementResults)[classIndex];

                std::list<size_t>::const_iterator deletedPointsIterator = std::begin(
                        refinementResult.deletedGridPointsIndexes);
                std::list<size_t>::const_iterator deletedPointsEnd = std::end(
                        refinementResult.deletedGridPointsIndexes);
                sendRefinementUpdates<size_t>(classIndex, DELETED_GRID_POINTS_LIST,
                                              deletedPointsIterator,
                                              deletedPointsEnd);
                std::list<sgpp::base::DataVector>::const_iterator addedPointsIterator = std::begin(
                        refinementResult.addedGridPoints);
                std::list<sgpp::base::DataVector>::const_iterator addedPointsEnd = std::end(
                        refinementResult.addedGridPoints);
                sendRefinementUpdates<sgpp::base::DataVector>(classIndex, ADDED_GRID_POINTS_LIST, addedPointsIterator,
                                                              addedPointsEnd);
            }
        }

        template<typename DataType>
        void MPIMethods::sendRefinementUpdates(size_t &classIndex, const RefinementResultsUpdateType updateType,
                                               typename std::list<DataType>::const_iterator &iterator,
                                               typename std::list<DataType>::const_iterator &listEnd) {
            while (iterator != listEnd) {
                    MPI_Packet *mpiPacket = new MPI_Packet;
                    RefinementResultNetworkMessage *networkMessage = (RefinementResultNetworkMessage *) mpiPacket->payload;

                    networkMessage->classIndex = classIndex;

                networkMessage->updateType = updateType;

                size_t numPointsInBuffer;
                switch (updateType) {
                    case DELETED_GRID_POINTS_LIST: {
                        numPointsInBuffer = fillBufferWithData<DataType>((void *) networkMessage->payload,
                                                                         (void *) std::end(networkMessage->payload),
                                                                         iterator,
                                                                         listEnd);
                    }
                        break;
                    case ADDED_GRID_POINTS_LIST:
                        numPointsInBuffer = fillBufferWithVectorData<DataType>((void *) networkMessage->payload,
                                                                               (void *) std::end(
                                                                                       networkMessage->payload),
                                                                               iterator,
                                                                               listEnd);
                        break;
                    }

                networkMessage->listLength = numPointsInBuffer;

                printf("Sending class %i update with %i deleted grid points", networkMessage->classIndex,
                       networkMessage->listLength);

                    PendingMPIRequest pendingMPIRequest;
                    pendingMPIRequest.buffer = mpiPacket;
                    pendingMPIRequests.push_back(pendingMPIRequest);

                //TODO: NiceToHave Send the smallest packet possible
                MPI_Ibcast(&networkMessage, sizeof(RefinementResultNetworkMessage), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                               MPI_COMM_WORLD, &(pendingMPIRequest.request));

                }
        }

        //TODO: This was imported from Merge
        size_t receiveMergeGridNetworkMessage(int gridversion, MergeGridNetworkMessage &networkMessage,
                                              base::DataVector &alphaVector) {
            if (gridversion != networkMessage.gridversion) {
                sgpp::base::application_exception applicationException(
                        "Received grid merge request with incorrect grid version!");
                throw applicationException;
            }

            double *payload = (double *) networkMessage.payload;
            for (size_t index; index < networkMessage.payloadLength; index++) {
                alphaVector[networkMessage.payloadOffset + index] = payload[index];
            }

            printf("Updated alpha values from network message offset %i, length %i", networkMessage.payloadOffset,
                   networkMessage.payloadLength);

        }

        //TODO: This was imported from Merge
/*
        size_t MPIMethods::sendRefinementResultPacket(size_t classIndex, RefinementResultsUpdateType updateType,
                                                      const RefinementResult &refinementResult, int offset,
                                                      std::list::iterator &iterator) {
            MPI_Packet *mpiPacket = new MPI_Packet;
            RefinementResultNetworkMessage *networkMessage = (RefinementResultNetworkMessage *) mpiPacket->payload;

            networkMessage->classIndex = classIndex;
            networkMessage->gridversion = refinementResult.gridversion;
            networkMessage->updateType = updateType;

            //Minimum between remaining grid points and maximum number of grid points that will still fit
            size_t numPointsInPacket = 0;

            if (updateType == DELETED_GRID_POINTS_LIST) {
                size_t endOfPayload = sizeof(RefinementResultNetworkMessage::payload) / sizeof(int);

                //Copy data from list into network buffer
                int *deletedGridPoints = (int *) networkMessage->payload;
                for (int bufferIndex = 0; bufferIndex < endOfPayload; bufferIndex++) {
                    deletedGridPoints[bufferIndex] = *iterator;
                    iterator++;
                    numPointsInPacket++;
                }
            } else if (updateType == ADDED_GRID_POINTS_LIST) {
                //TODO: These vectors do not have a constant grid size
                size_t endOfPayload = sizeof(RefinementResultNetworkMessage::payload) / sizeof(sgpp::base::DataVector);
                //Copy data from list into network buffer
                sgpp::base::DataVector *addedGridPoints = (sgpp::base::DataVector *) networkMessage->payload;
                for (int bufferIndex = 0; bufferIndex < endOfPayload; bufferIndex++) {
                    addedGridPoints[bufferIndex] = *iterator;
                    iterator++;
                    numPointsInPacket++;
                }
            }

            networkMessage->listLength = numPointsInPacket;

            printf("Sending update for class %zu with %i modifications", classIndex,
                   networkMessage->listLength);
            PendingMPIRequest pendingMPIRequest;
            pendingMPIRequest.buffer = mpiPacket;
            pendingMPIRequests.push_back(pendingMPIRequest);

            //Send the smallest packet possible
            MPI_Ibcast(&networkMessage, numPointsInPacket + 3, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                       MPI_COMM_WORLD, &(pendingMPIRequest.request));
            return numPointsInPacket;
        }
*/

        //TODO: This was imported from Merge
        size_t
        sendMergeGridNetworkMessage(size_t classIndex, base::DataVector &alphaVector, size_t offset, size_t length) {
            //TODO
            MPI_Packet *mpiPacket = new MPI_Packet;
            MergeGridNetworkMessage *networkMessage = (MergeGridNetworkMessage *) mpiPacket->payload;

            networkMessage->classIndex = classIndex;
//            networkMessage->gridversion = gridversion;
            networkMessage->payloadOffset = offset;

            size_t endOfPayload = sizeof(MergeGridNetworkMessage::payload) / sizeof(double);
            size_t numPointsInPacket = 0;

            size_t endOfCopy = std::min(length, endOfPayload);

            double *alphaValues = (double *) networkMessage->payload;
            for (size_t bufferIndex = 0; bufferIndex < endOfCopy; bufferIndex++) {
                alphaValues[bufferIndex] = alphaVector[offset + bufferIndex];
                numPointsInPacket++;
            }

            printf("Sending merge for class %zu offset %zu and %zu modifications", classIndex, offset,
                   numPointsInPacket);
            networkMessage->payloadLength = numPointsInPacket;
            return numPointsInPacket;
        }


        //TODO: !!!!!!!!! Compiler Errors

        //TODO: Ensure compiler calls the correct method
        template<typename DataType>
        size_t
        MPIMethods::fillBufferWithVectorData(void *buffer, void *bufferEnd,
                                             typename std::vector<DataType>::const_iterator iterator,
                                             typename DataType::const_iterator listEnd) {
            //TODO: Implement vector
            DataType *bufferPointer = (DataType *) buffer;
            size_t copiedVectors = 0;
            while (bufferPointer + iterator->size() * sizeof(DataType) < bufferEnd && iterator != listEnd) {
                std::memcpy(&((*iterator)[0]), bufferPointer, iterator->size());
                iterator++;
                copiedVectors++;
            }
            return copiedVectors;
        }

        template<typename DataType>
        size_t
        MPIMethods::fillBufferWithData(void *buffer, void *bufferEnd, typename DataType::const_iterator iterator,
                                       typename DataType::const_iterator listEnd) {
            DataType *bufferPointer = (DataType *) buffer;
            size_t copiedValues = 0;
            while (bufferPointer + sizeof(DataType) < bufferEnd && iterator != listEnd) {
                *bufferPointer = *iterator;
                iterator++;
                copiedValues++;
            }
            return copiedValues;
        }

        void MPIMethods::processCompletedMPIRequests() {
            MPI_Status mpiStatus;
            int operationCompleted;

            for (std::vector<sgpp::datadriven::PendingMPIRequest>::iterator pendingMPIRequestIterator = pendingMPIRequests.begin();
                 pendingMPIRequestIterator != pendingMPIRequests.end();
                 pendingMPIRequestIterator++) {
                MPI_Test(&(pendingMPIRequestIterator->request), &operationCompleted, &mpiStatus);
                if (operationCompleted) {
                    //Execute the callback
                    pendingMPIRequestIterator->callback(*pendingMPIRequestIterator);

                    if (pendingMPIRequestIterator->disposeAfterCallback) {
                        //TODO Deleting a void pointer here
                        delete[] pendingMPIRequestIterator->buffer;
                        pendingMPIRequests.erase(pendingMPIRequestIterator);
                    }
                }
            }
        }

        void MPIMethods::waitForMPIRequestsToComplete() {
            for (PendingMPIRequest &pendingMPIRequest : pendingMPIRequests) {
                MPI_Wait(&(pendingMPIRequest.request), MPI_STATUS_IGNORE);
            }

            processCompletedMPIRequests();
        }

        void MPIMethods::receiveGridComponentsUpdate(LearnerSGDEOnOffParallel *learnerInstance,
                                                     RefinementResultNetworkMessage *networkMessage) {
            //TODO
            RefinementResult refinementResult;

            void *bufferEnd = std::end(networkMessage->payload);
            switch (networkMessage->updateType) {
                case DELETED_GRID_POINTS_LIST: {
                    //TODO: This probably casts each unsigned char into a size_t instead of re-interpreting
                    size_t *bufferIterator = (size_t *) networkMessage->payload;
                    while (bufferIterator < bufferEnd) {
                        refinementResult.deletedGridPointsIndexes.push_back(*bufferIterator);
                        bufferIterator++;
                    }
                }
                    break;
                case ADDED_GRID_POINTS_LIST: {
                    //TODO double?
                    double *bufferIterator = (double *) networkMessage->payload;
                    size_t dimensionality = learnerInstance->getDimensionality();
                    while (bufferIterator < bufferEnd) {
                        size_t index = 0;
                        sgpp::base::DataVector dataVector;
                        while (index < dimensionality && bufferIterator < bufferEnd) {
                            dataVector.push_back(*bufferIterator);
                            index++;
                            bufferIterator++;
                        }
                        refinementResult.addedGridPoints.push_back(dataVector);
                    }
                    break;
                }
            }
            learnerInstance->updateVariablesAfterRefinement(&refinementResult, networkMessage->classIndex,
                                                            learnerInstance->getDensityFunctions()[networkMessage->classIndex].first.get());
        }

        void MPIMethods::processIncomingMPICommands(LearnerSGDEOnOffParallel *learnerInstance, MPI_Packet *mpiPacket) {
            switch (mpiPacket->commandID) {
                case START_SYNCHRONIZE_PACKETS:
                    startSynchronizingPackets();
                    break;
                case UPDATE_GRID:
                    receiveGridComponentsUpdate(learnerInstance,
                                                (RefinementResultNetworkMessage *) (mpiPacket->payload));
                    break;
                case MERGE_GRID:

                    break;
                case ASSIGN_BATCH:
                    break;
                default:
                    printf("Error: MPI unknown command id: %i", mpiPacket->commandID);
                    break;
            }
        }

        void MPIMethods::startSynchronizingPackets() {
            //TODO
        }

    }
}