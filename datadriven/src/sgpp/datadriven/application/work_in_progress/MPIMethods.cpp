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
#include <thread>

namespace sgpp {
    namespace datadriven {

        //Pending MPI Requests
        std::list<sgpp::datadriven::PendingMPIRequest> MPIMethods::pendingMPIRequests;
        int MPIMethods::mpiWorldSize = -1;
        LearnerSGDEOnOffParallel *MPIMethods::learnerInstance;


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
                auto *mpiPacket = new MPI_Packet;
                auto &unicastInputRequest = createPendingMPIRequest(mpiPacket);
                unicastInputRequest.disposeAfterCallback = false;
                unicastInputRequest.callback = [](PendingMPIRequest &request) {
                    std::cout << "Incoming MPI unicast" << std::endl;
                    processIncomingMPICommands(request.buffer);

//                    std::cout << "Zeroing MPI Request" << std::endl;
//                    std::memset(request.request, 0, sizeof(MPI_Request));

                    MPI_Irecv(request.buffer, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
                              MPI_ANY_TAG, MPI_COMM_WORLD, &(request.request));
                };

                MPI_Irecv(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                          MPI_COMM_WORLD,
                          &(unicastInputRequest.request));

                std::cout << "Started listening for unicasts from any sources" << std::endl;
            }
            if (!isMaster()) {
                auto *mpiPacket = new MPI_Packet;
                auto &broadcastInputRequest = createPendingMPIRequest(mpiPacket);
                broadcastInputRequest.disposeAfterCallback = false;
                broadcastInputRequest.callback = [](PendingMPIRequest &request) {
                    std::cout << "Incoming MPI broadcast" << std::endl;
                    processIncomingMPICommands(request.buffer);

//                    std::cout << "Zeroing MPI Request" << std::endl;
//                    std::memset(request.request, 0, sizeof(MPI_Request));

                    MPI_Ibcast(request.buffer, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                               MPI_COMM_WORLD, &(request.request));
                };
                MPI_Ibcast(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK, MPI_COMM_WORLD,
                           &(broadcastInputRequest.request));

                MPIMethods::learnerInstance = learnerInstance;
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
                          << " with " << networkMessage->listLength
                          << " modifications" << std::endl;

                sendIBcast(mpiPacket);
            }
        }

        void MPIMethods::bcastCommandNoArgs(MPI_COMMAND_ID commandId) {
            MPI_Packet *mpiPacket = new MPI_Packet;
            mpiPacket->commandID = commandId;

            sendIBcast(mpiPacket);
        }

        void MPIMethods::sendCommandNoArgs(const int destinationRank, MPI_COMMAND_ID commandId) {
            MPI_Packet *mpiPacket = new MPI_Packet;
            mpiPacket->commandID = commandId;

            sendISend(destinationRank, mpiPacket);
        }

        void MPIMethods::sendIBcast(MPI_Packet *mpiPacket) {
            PendingMPIRequest &pendingMPIRequest = createPendingMPIRequest(mpiPacket);

            MPI_Ibcast(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                       MPI_COMM_WORLD, &(pendingMPIRequest.request));
        }

        void MPIMethods::sendISend(const int destinationRank, MPI_Packet *mpiPacket) {
            PendingMPIRequest &pendingMPIRequest = createPendingMPIRequest(mpiPacket);

            //Point to the request in vector instead of stack
            MPI_Isend(mpiPacket, sizeof(MPI_Packet), MPI_UNSIGNED_CHAR, destinationRank, COMMAND_TAG,
                      MPI_COMM_WORLD, &(pendingMPIRequest.request));
        }

        PendingMPIRequest &MPIMethods::createPendingMPIRequest(MPI_Packet *mpiPacket) {
            pendingMPIRequests.emplace_back();
            PendingMPIRequest &pendingMPIRequest = pendingMPIRequests.back();
            pendingMPIRequest.disposeAfterCallback = true;
            pendingMPIRequest.callback = [](PendingMPIRequest &request) {
                std::cout << "Pending MPI request " << &request << "completed. " << std::endl;
            };
            pendingMPIRequest.buffer = mpiPacket;
            return pendingMPIRequest;
        }

        //TODO: This was imported from Merge
        size_t MPIMethods::receiveMergeGridNetworkMessage(int gridversion, MergeGridNetworkMessage &networkMessage,
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
            //TODO: DANGEROUS
            pendingMPIRequests.push_back(pendingMPIRequest);

            //Send the smallest packet possible
            MPI_Ibcast(&networkMessage, numPointsInPacket + 3, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                       MPI_COMM_WORLD, &(pendingMPIRequest.request));
            return numPointsInPacket;
        }
*/

        size_t
        MPIMethods::sendMergeGridNetworkMessage(size_t classIndex, base::DataVector &alphaVector) {
            size_t offset = 0;
            while (offset < alphaVector.size()) {
                MPI_Packet *mpiPacket = new MPI_Packet;
                mpiPacket->commandID = MERGE_GRID;

                MergeGridNetworkMessage *networkMessage = (MergeGridNetworkMessage *) mpiPacket->payload;

                networkMessage->classIndex = classIndex;
//            networkMessage->gridversion = gridversion;
                networkMessage->payloadOffset = offset;

                size_t endOfPayload = sizeof(MergeGridNetworkMessage::payload) / sizeof(double);
                size_t numPointsInPacket = 0;

                size_t endOfCopy = std::min(alphaVector.size(), endOfPayload);

                double *alphaValues = (double *) networkMessage->payload;
                for (size_t bufferIndex = 0; bufferIndex < endOfCopy; bufferIndex++) {
                    alphaValues[bufferIndex] = alphaVector[offset + bufferIndex];
                    numPointsInPacket++;
                }
                networkMessage->payloadLength = numPointsInPacket;

                std::cout << "Sending merge for class " << classIndex
                          << " offset " << offset
                          << " and " << numPointsInPacket << " values" << std::endl;
                sendISend(MPI_MASTER_RANK, mpiPacket);
                offset += numPointsInPacket;
            }

            return offset;
        }


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

            std::cout << "Checking " << pendingMPIRequests.size() << " pending MPI requests" << std::endl;

            for (auto pendingMPIRequestIterator = pendingMPIRequests.begin();
                 pendingMPIRequestIterator != pendingMPIRequests.end();
                 pendingMPIRequestIterator++) {
                MPI_Status mpiStatus{};
                int operationCompleted;

                std::cout << "Testing request " << &*pendingMPIRequestIterator << std::endl;
                if (MPI_Test(&(pendingMPIRequestIterator->request), &operationCompleted, &mpiStatus)) {
                    std::cout << "Error MPI Test reported" << std::endl;
                    exit(-1);
                }
                std::cout << "Pending request has status " << operationCompleted
                          << " MPI_ERROR: " << mpiStatus.MPI_ERROR
                          << " MPI SOURCE: " << mpiStatus.MPI_SOURCE
                          << " MPI TAG: " << mpiStatus.MPI_TAG << std::endl;
                if (operationCompleted != 0) {

                    std::cout << "Executing callback" << std::endl;
                    //Execute the callback
                    pendingMPIRequestIterator->callback(*pendingMPIRequestIterator);
                    std::cout << "Callback complete" << std::endl;


                    if (pendingMPIRequestIterator->disposeAfterCallback) {
                        //TODO Deleting a void pointer here
                        std::cout << "Attempting to delete pending mpi request" << std::endl;
                        delete[] pendingMPIRequestIterator->buffer;
                        //TODO: !!! DELETING STUFF HERE WILL MOVE OTHER PENDING REQUESTS. THIS CAUSES SEGMENTATION FAULTS
                        pendingMPIRequests.erase(pendingMPIRequestIterator);
                        std::cout << "Deleted pending mpi request" << std::endl;
                    } else {
//                        std::cout << "Zeroing MPI Request" << std::endl;
//                        std::memset(pendingMPIRequestIterator->request, 0, sizeof(MPI_Request));
//
//                        std::cout << "Zeroing Buffer" << std::endl;
//                        std::memset(pendingMPIRequestIterator->buffer, 0, sizeof(MPI_Packet));

                    }
                    std::cout << "Relaunching processCompletedMPIRequests" << std::endl;
                    processCompletedMPIRequests();
                    break;
                }
            }
        }

        void MPIMethods::waitForAllMPIRequestsToComplete() {
            for (PendingMPIRequest &pendingMPIRequest : pendingMPIRequests) {
                MPI_Wait(&(pendingMPIRequest.request), MPI_STATUS_IGNORE);
            }
        }

        void MPIMethods::waitForAnyMPIRequestsToComplete() {
            int completedRequest;
            std::cout << "Not implemented" << std::endl;
            exit(-1);
            //TODO: The pending MPI requests aren't held in memory sequentially
//            MPI_Waitany(pendingMPIRequests.size(), &pendingMPIRequests[0].request, &completedRequest, MPI_STATUS_IGNORE);
        }

        void MPIMethods::receiveGridComponentsUpdate(RefinementResultNetworkMessage *networkMessage) {
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

        void MPIMethods::processIncomingMPICommands(MPI_Packet *mpiPacket) {
            std::cout << "Processing incoming command " << mpiPacket->commandID << std::endl;
            switch (mpiPacket->commandID) {
                case START_SYNCHRONIZE_PACKETS:
                    startSynchronizingPackets();
                    break;
                case END_SYNCHRONIZE_PACKETS:
                    endSynchronizingPackets();
                    break;
                case UPDATE_GRID:
                    receiveGridComponentsUpdate((RefinementResultNetworkMessage *) (mpiPacket->payload));
                    break;
                case MERGE_GRID:
                    std::cout << "Merge grid not implemented" << std::endl;
//                    receiveMergeGridNetworkMessage(learnerInstance->getGridVersion(), reinterpret_cast<MergeGridNetworkMessage>(mpiPacket->payload), learnerInstance->getAlphaVector());
                    break;
                case ASSIGN_BATCH:
                    runBatch(mpiPacket);
                    break;
                case SHUTDOWN:
                    std::cout << "Worker shutdown requested" << std::endl;
                    learnerInstance->shutdown();
                    break;
                case NULL_COMMAND:
                    std::cout << "Error: Incoming command has undefined command id" << std::endl;
                    exit(-1);
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

        void MPIMethods::assignBatch(const int workerID, size_t batchOffset, size_t batchSize, bool doCrossValidation) {
            MPI_Packet *mpiPacket = new MPI_Packet;
            mpiPacket->commandID = ASSIGN_BATCH;

            auto *message = (AssignBatchNetworkMessage *) mpiPacket->payload;
            message->batchOffset = batchOffset;
            message->batchSize = batchSize;
            message->doCrossValidation = doCrossValidation;

            sendISend(workerID, mpiPacket);
        }

        void MPIMethods::runBatch(MPI_Packet *mpiPacket) {
            auto *message = (AssignBatchNetworkMessage *) mpiPacket->payload;
            std::cout << "runbatch dim " << learnerInstance->getDimensionality() << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(1));
            std::cout << "creating dataset" << std::endl;
            Dataset dataset{message->batchSize, learnerInstance->getDimensionality()};
            learnerInstance->workBatch(dataset, message->batchOffset, message->doCrossValidation);
        }

        int MPIMethods::getWorldSize() {
            return mpiWorldSize;
        }

    }
}