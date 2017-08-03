//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_NETWORKMESSAGEDATA_H
#define SGPP_NETWORKMESSAGEDATA_H

#include <mpi.h>
#include <functional>

namespace sgpp {
    namespace datadriven {
        static const int MPI_PACKET_MAX_PAYLOAD_SIZE = 512;
        static const int MPI_MASTER_RANK = 0;
        static const int MPI_MAX_PROCESSOR_NAME_LENGTH = 256;

        enum MPI_COMMAND_ID {
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
            unsigned int gridversion;
            unsigned int classIndex;
            RefinementResultsUpdateType updateType;
            unsigned int listLength;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE - 2 * sizeof(unsigned int) -
                                   sizeof(RefinementResultsUpdateType))];
        };

        struct MergeGridNetworkMessage {
            unsigned int gridversion;
            unsigned int classIndex;
            unsigned int payloadOffset;
            unsigned int payloadLength;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE - 4 * sizeof(int))];
        };

        struct AssignBatchNetworkMessage {
            unsigned int batchOffset;
            unsigned int batchSize;
            bool doCrossValidation;
        };

    }
}

#endif //SGPP_NETWORKMESSAGEDATA_H
