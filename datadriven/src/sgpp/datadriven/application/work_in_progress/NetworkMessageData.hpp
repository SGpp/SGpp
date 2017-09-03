//
// Created by Vincent_Bode on 29.06.2017.
//

#ifndef SGPP_NETWORKMESSAGEDATA_H
#define SGPP_NETWORKMESSAGEDATA_H

#define MPI_PACKET_MAX_PAYLOAD_SIZE 4096
#define MPI_MASTER_RANK 0
#define MPI_MAX_PROCESSOR_NAME_LENGTH 256
#define MPI_TAG_HIGH_PRIORITY_NO_BLOCK 42
#define MPI_TAG_STANDARD_COMMAND 41

#define REFINENEMT_RESULT_PAYLOAD_SIZE (MPI_PACKET_MAX_PAYLOAD_SIZE\
                                   - 3 * sizeof(unsigned long)\
                                   - sizeof(RefinementResultsUpdateType))

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
            UPDATE_SYSTEM_MATRIX_DECOMPOSITION,
            SHUTDOWN,
            WORKER_SHUTDOWN_SUCCESS
        };

        struct MPI_Packet {
            MPI_COMMAND_ID commandID;
            unsigned char payload[MPI_PACKET_MAX_PAYLOAD_SIZE];
        };

        enum RefinementResultsUpdateType {
            ADDED_GRID_POINTS_LIST,
            DELETED_GRID_POINTS_LIST,
            CHOLESKY_DECOMPOSITION
        };

        struct RefinementResultNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            unsigned long listLength;
            RefinementResultsUpdateType updateType;

            unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE];
        };

        struct RefinementResultCholeskyNetworkMessage {
            unsigned long matrixWidth;
            unsigned long matrixHeight;
            unsigned long offset;

            unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE - 3 * sizeof(unsigned long)];
        };

        struct MergeGridNetworkMessage {
            unsigned long gridversion;
            unsigned long classIndex;
            unsigned long payloadOffset;
            unsigned long payloadLength;
            unsigned long batchSize;
            unsigned long batchOffset;
            unsigned long alphaTotalSize;

            unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
                                   - 7 * sizeof(unsigned long))];
        };

        struct AssignBatchNetworkMessage {
            unsigned long batchOffset;
            unsigned long batchSize;
            bool doCrossValidation;
        };

        struct AssignSystemMatrixUpdateNetworkMessage {
            unsigned long classIndex;
            unsigned long gridversion;
        };

        static_assert(sizeof(size_t) <= sizeof(unsigned long), "size_t larger than unsigned long");
        static_assert(sizeof(MergeGridNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
                      "Merge Grid Network Message too long.");
        static_assert(sizeof(RefinementResultNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
                      "Refinement result Network Message too long.");
        static_assert(sizeof(RefinementResultCholeskyNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
                      "Refinement result Cholesky Network Message too long.");


    }
}

#endif //SGPP_NETWORKMESSAGEDATA_H
