//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_NETWORKMESSAGEDATA_H
#define SGPP_NETWORKMESSAGEDATA_H

#define MPI_PACKET_MAX_PAYLOAD_SIZE 512
#define MPI_MASTER_RANK 0
#define MPI_MAX_PROCESSOR_NAME_LENGTH 256

#include <mpi.h>
#include <functional>

namespace sgpp {
    namespace datadriven {

        enum MPI_COMMAND_ID {
            NULL_COMMAND,
            START_SYNCHRONIZE_PACKETS,
            END_SYNCHRONIZE_PACKETS,
            UPDATE_GRID,
            MERGE_GRID,
            ASSIGN_BATCH,
            SHUTDOWN
        };

        enum MPI_COMMAND_TAG {
            COMMAND_TAG,
        };

        struct MPI_Packet {
            MPI_COMMAND_ID commandID;
            unsigned char payload[MPI_PACKET_MAX_PAYLOAD_SIZE];
        };

        struct PendingMPIRequest {
            MPI_Request request;
            sgpp::datadriven::MPI_Packet *buffer;
            std::function<void(PendingMPIRequest &)> callback;
            bool disposeAfterCallback;
        };

        enum RefinementResultsUpdateType {
            ADDED_GRID_POINTS_LIST,
            DELETED_GRID_POINTS_LIST
        };

        struct RefinementResultNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            RefinementResultsUpdateType updateType;
            unsigned long listLength;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
                                   - 3 * sizeof(unsigned long)
                                   - sizeof(RefinementResultsUpdateType))];
        };

        struct MergeGridNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            unsigned long payloadOffset;
            unsigned long payloadLength;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
                                   - 4 * sizeof(unsigned long))];
        };

        struct AssignBatchNetworkMessage {
            unsigned long batchOffset;
            unsigned long batchSize;
            bool doCrossValidation;
        };

    }
}

#endif //SGPP_NETWORKMESSAGEDATA_H
