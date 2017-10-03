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
            ASSIGN_BATCH
        };

        enum MPI_COMMAND_TAG {

        };

        struct PendingMPIRequest {
            MPI_Request request;
            MPI_Packet *buffer;
            std::function<void(void *)> callback;
            bool disposeAfterCallback;
        };

        struct MPI_Packet {
            MPI_COMMAND_ID commandID;
            unsigned char payload[MPI_PACKET_MAX_PAYLOAD_SIZE];
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
            int gridversion;
            int classIndex;
            int payloadOffset;
            int payloadLength;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE - 4 * sizeof(int))];
        };

    }
}

#endif //SGPP_NETWORKMESSAGEDATA_H
