//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_NETWORKMESSAGEDATA_H
#define SGPP_NETWORKMESSAGEDATA_H

#define MPI_PACKET_MAX_PAYLOAD_SIZE 4096
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

        enum RefinementResultsUpdateType {
            ADDED_GRID_POINTS_LIST,
            DELETED_GRID_POINTS_LIST
        };

        struct RefinementResultNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            unsigned long listLength;
            RefinementResultsUpdateType updateType;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
                                   - 3 * sizeof(unsigned long)
                                   - sizeof(RefinementResultsUpdateType))];
        };

        struct MergeGridNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            unsigned long payloadOffset;
            unsigned long payloadLength;
            unsigned long batchSize;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
                                   - 5 * sizeof(unsigned long))];
        };

        struct AssignBatchNetworkMessage {
            unsigned long batchOffset;
            unsigned long batchSize;
            bool doCrossValidation;
        };

        static_assert(sizeof(size_t) <= sizeof(unsigned long), "size_t larger than unsigned long");
        static_assert(sizeof(MergeGridNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
                      "Merge Grid Network Message too long.");
        static_assert(sizeof(RefinementResultNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
                      "Refinement result Network Message too long.");


    }
}

#endif //SGPP_NETWORKMESSAGEDATA_H
