//
// Created by Vincent_Bode on 29.06.2017.
//

#include <ctime>
#include <mpi.h>
#include <cstring>

#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
    namespace datadriven {

        //Pending MPI Requests
        std::vector<sgpp::datadriven::PendingMPIRequest> MPIMethods::pendingMPIRequests;
        int MPIMethods::mpiWorldSize = -1;


        bool MPIMethods::isMaster() {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            return rank == MPI_MASTER_RANK;
        }

        void MPIMethods::initMPI(LearnerSGDEOnOffParallel *learnerInstance) {
            MPI_Init(nullptr, nullptr);

            // Get World Size
            MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

            //Get Rank
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            //Get Processor Name
            char mpiProcessorName[MPI_MAX_PROCESSOR_NAME_LENGTH];
            int nameLength;
            MPI_Get_processor_name(mpiProcessorName, &nameLength);

            std::cout << "Processor " << mpiProcessorName << " (rank " << world_rank
                      << ") has joined MPI pool of size " << mpiWorldSize << std::endl;

            //Setup receiving messages from master/workers
            {
                sgpp::datadriven::PendingMPIRequest unicastInputRequest;
                auto *mpiPacket = new MPI_Packet;
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
                std::cout << "Started listening for unicasts from any sources" << std::endl;
            }
            if (!isMaster()) {
                PendingMPIRequest broadcastInputRequest;
                auto *mpiPacket = new MPI_Packet;
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

                std::cout << "Started listening for broadcasts from task master" << std::endl;
            }

        }

        void MPIMethods::synchronizeBarrier() {
            MPI_Barrier(MPI_COMM_WORLD);
        }

        void MPIMethods::sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults) {

            std::cout << "Updating the grid components on workers..." << std::endl;

            for (size_t classIndex = 0; classIndex < refinementResults->size(); classIndex++) {
                RefinementResult refinementResult = (*refinementResults)[classIndex];

                std::cout << "Updating grid for class " << classIndex
                          << " (" << refinementResult.addedGridPoints.size() << " additions, "
                          << refinementResult.deletedGridPointsIndexes.size() << " deletions)" << std::endl;

                std::list<size_t>::const_iterator deletedPointsIterator = std::begin(
                        refinementResult.deletedGridPointsIndexes);
                std::list<size_t>::const_iterator deletedPointsEnd = std::end(
                        refinementResult.deletedGridPointsIndexes);
                sendRefinementUpdates<std::list<size_t>::const_iterator>(classIndex, DELETED_GRID_POINTS_LIST,
                                                                         deletedPointsIterator,
                                                                         deletedPointsEnd);
                std::list<sgpp::base::DataVector>::const_iterator addedPointsIterator = std::begin(
                        refinementResult.addedGridPoints);
                std::list<sgpp::base::DataVector>::const_iterator addedPointsEnd = std::end(
                        refinementResult.addedGridPoints);
                sendRefinementUpdates<std::list<DataVector>::const_iterator>(classIndex, ADDED_GRID_POINTS_LIST,
                                                                             addedPointsIterator,
                                                                             addedPointsEnd);
            }

            std::cout << "Finished updating the grid components on workers." << std::endl;
        }

        template<typename Iterator>
        void MPIMethods::sendRefinementUpdates(size_t &classIndex, const RefinementResultsUpdateType updateType,
                                               Iterator &iterator,
                                               Iterator &listEnd) {
            while (iterator != listEnd) {
                auto *mpiPacket = new MPI_Packet;
                auto *networkMessage = (RefinementResultNetworkMessage *) mpiPacket->payload;

                    networkMessage->classIndex = classIndex;

                networkMessage->updateType = updateType;

                size_t numPointsInBuffer;
                switch (updateType) {
                    case DELETED_GRID_POINTS_LIST: {
                        numPointsInBuffer = fillBufferWithData<Iterator>((void *) networkMessage->payload,
                                                                         (void *) std::end(networkMessage->payload),
                                                                         iterator,
                                                                         listEnd);
                    }
                        break;
                    case ADDED_GRID_POINTS_LIST:
                        numPointsInBuffer = fillBufferWithVectorData<Iterator, double_t>(
                                (void *) networkMessage->payload,
                                (void *) std::end(
                                                                                       networkMessage->payload),
                                iterator,
                                listEnd,
                                sizeof(double)); //TODO: Size of double represents the size of a data value in DataVector
                        break;
                    default:
                        std::cout << "ERROR: Unknown update type" << std::endl;
                        exit(-1);
                    }

                networkMessage->listLength = numPointsInBuffer;

                std::cout << "Sending updated for class " << networkMessage->classIndex
                          << "with " << networkMessage->listLength
                          << " modifications" << std::endl;

                sendIBcast(mpiPacket);
            }
        }

        void MPIMethods::bcastCommandNoArgs(MPI_COMMAND_ID commandId) {
            MPI_Packet mpiPacket{};
            mpiPacket.commandID = commandId;

            sendIBcast(&mpiPacket);
        }

        void MPIMethods::sendCommandNoArgs(int destinationRank, MPI_COMMAND_ID commandId) {
            MPI_Packet mpiPacket{};
            mpiPacket.commandID = commandId;

            sendISend(destinationRank, &mpiPacket);
        }

        void MPIMethods::sendIBcast(MPI_Packet *mpiPacket) {
            PendingMPIRequest pendingMPIRequest;
            pendingMPIRequest.buffer = mpiPacket;
            pendingMPIRequests.push_back(pendingMPIRequest);

            MPI_Ibcast(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                       MPI_COMM_WORLD, &(pendingMPIRequest.request));
        }

        void MPIMethods::sendISend(int destinationRank, MPI_Packet *mpiPacket) {
            PendingMPIRequest pendingMPIRequest;
            pendingMPIRequest.buffer = mpiPacket;
            pendingMPIRequests.push_back(pendingMPIRequest);

            MPI_Isend(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, destinationRank, MPI_ANY_TAG,
                      MPI_COMM_WORLD, &(pendingMPIRequest.request));

        }

        //TODO: This was imported from Merge
        size_t receiveMergeGridNetworkMessage(int gridversion, MergeGridNetworkMessage &networkMessage,
                                              base::DataVector &alphaVector) {
            if (gridversion != networkMessage.gridversion) {
                sgpp::base::application_exception applicationException(
                        "Received grid merge request with incorrect grid version!");
                throw applicationException;
            }

            auto *payload = (double *) networkMessage.payload;
            for (size_t index = 0; index < networkMessage.payloadLength; index++) {
                alphaVector[networkMessage.payloadOffset + index] = payload[index];
            }

            std::cout << "Updated alpha values from network message offset " << networkMessage.payloadOffset
                      << ", length " << networkMessage.payloadLength;
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

//        //TODO: This was imported from Merge
//        size_t
//        sendMergeGridNetworkMessage(size_t classIndex, base::DataVector &alphaVector, size_t offset, size_t length) {
//            //TODO
//            MPI_Packet *mpiPacket = new MPI_Packet;
//            MergeGridNetworkMessage *networkMessage = (MergeGridNetworkMessage *) mpiPacket->payload;
//
//            networkMessage->classIndex = classIndex;
////            networkMessage->gridversion = gridversion;
//            networkMessage->payloadOffset = offset;
//
//            size_t endOfPayload = sizeof(MergeGridNetworkMessage::payload) / sizeof(double);
//            size_t numPointsInPacket = 0;
//
//            size_t endOfCopy = std::min(length, endOfPayload);
//
//            double *alphaValues = (double *) networkMessage->payload;
//            for (size_t bufferIndex = 0; bufferIndex < endOfCopy; bufferIndex++) {
//                alphaValues[bufferIndex] = alphaVector[offset + bufferIndex];
//                numPointsInPacket++;
//            }
//
//            printf("Sending merge for class %zu offset %zu and %zu modifications", classIndex, offset,
//                   numPointsInPacket);
//            networkMessage->payloadLength = numPointsInPacket;
//            return numPointsInPacket;
//        }


        //TODO: !!!!!!!!! Compiler Errors

        //TODO: Ensure compiler calls the correct method
        template<typename Iterator, typename ValueType>
        size_t
        MPIMethods::fillBufferWithVectorData(void *buffer, const void *bufferEnd,
                                             Iterator &iterator,
                                             Iterator &listEnd, size_t sizeOfDataType) {
            //TODO: Implement vector
            auto *bufferPointer = buffer;
            size_t copiedVectors = 0;
            while (iterator != listEnd) {
                auto dataVector = (typename std::vector<ValueType>) *iterator;
                size_t vectorMemLength = dataVector.size() * sizeOfDataType;

                if (bufferPointer + vectorMemLength >= bufferEnd) {
                    break;
                }

                std::memcpy(&(dataVector[0]), bufferPointer, vectorMemLength);
                bufferPointer += vectorMemLength;
                iterator++;
                copiedVectors++;
            }
            return copiedVectors;
        }

        template<typename Iterator>
        size_t
        MPIMethods::fillBufferWithData(void *buffer, void *bufferEnd, Iterator &iterator,
                                       Iterator &listEnd) {
            auto *bufferPointer = (typename std::iterator_traits<Iterator>::value_type *) buffer;
            size_t copiedValues = 0;
            while (bufferPointer + sizeof(typename std::iterator_traits<Iterator>::value_type) < bufferEnd &&
                   iterator != listEnd) {
                *bufferPointer = *iterator;

                bufferPointer++;
                iterator++;
                copiedValues++;
            }
            return copiedValues;
        }

        void MPIMethods::processCompletedMPIRequests() {
            MPI_Status mpiStatus{};
            int operationCompleted;

            for (auto pendingMPIRequestIterator = pendingMPIRequests.begin();
                 pendingMPIRequestIterator != pendingMPIRequests.end();
                 pendingMPIRequestIterator++) {
                MPI_Test(&(pendingMPIRequestIterator->request), &operationCompleted, &mpiStatus);
                if (operationCompleted != 0) {
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
            RefinementResult refinementResult = sgpp::datadriven::RefinementResult();

            void *bufferEnd = std::end(networkMessage->payload);
            switch (networkMessage->updateType) {
                case DELETED_GRID_POINTS_LIST: {
                    //TODO: This probably casts each unsigned char into a size_t instead of re-interpreting
                    auto *bufferIterator = (size_t *) networkMessage->payload;
                    while (bufferIterator < bufferEnd) {
                        refinementResult.deletedGridPointsIndexes.push_back(*bufferIterator);
                        bufferIterator++;
                    }
                }
                    break;
                case ADDED_GRID_POINTS_LIST: {
                    //TODO double?
                    auto *bufferIterator = (double *) networkMessage->payload;
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
                case END_SYNCHRONIZE_PACKETS:
                    endSynchronizingPackets();
                    break;
                case UPDATE_GRID:
                    receiveGridComponentsUpdate(learnerInstance,
                                                (RefinementResultNetworkMessage *) (mpiPacket->payload));
                    break;
                case MERGE_GRID:
                    std::cout << "Merge grid not implemented" << std::endl;
                    break;
                case ASSIGN_BATCH:
                    std::cout << "Assign batch not implemented" << std::endl;
                    assignBatch(mpiPacket,
                                learnerInstance);
                    break;
                case SHUTDOWN:
                    std::cout << "Worker shutdown requested" << std::endl;
                    learnerInstance->shutdown();
                    break;
                default:
                    std::cout << "Error: MPI unknown command id: " << mpiPacket->commandID << std::endl;
                    exit(-1);
            }
        }

        void MPIMethods::startSynchronizingPackets() {
            //TODO
            std::cout << "Synchronizing not implemented" << std::endl;
        }

        void MPIMethods::finalizeMPI() {
            MPI_Finalize();
        }

        void MPIMethods::endSynchronizingPackets() {
            //TODO
            std::cout << "Synchronizing not implemented" << std::endl;
        }

        void MPIMethods::assignBatch(MPI_Packet *mpiPacket, LearnerSGDEOnOffParallel *learnerInstance) {
            auto *message = (AssignBatchNetworkMessage *) mpiPacket->payload;
            Dataset dataset{message->batchSize, learnerInstance->getDimensionality()};
            learnerInstance->assembleNextBatchData(&dataset, &(static_cast<size_t > (message->batchOffset)));
            //TODO continue
            learnerInstance->workBatch(dataset);
        }

    }
}