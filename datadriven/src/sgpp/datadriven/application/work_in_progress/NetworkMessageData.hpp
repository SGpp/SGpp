// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_NETWORKMESSAGEDATA_H
#define SGPP_NETWORKMESSAGEDATA_H

#define MPI_PACKET_MAX_PAYLOAD_SIZE 4096
#define MPI_MASTER_RANK 0
#define MPI_MAX_PROCESSOR_NAME_LENGTH 256
#define MPI_TAG_HIGH_PRIORITY_NO_BLOCK 42
#define MPI_TAG_STANDARD_COMMAND 41

#define REFINENEMT_RESULT_PAYLOAD_SIZE (MPI_PACKET_MAX_PAYLOAD_SIZE\
                                   - 3 * sizeof(size_t)\
                                   - sizeof(RefinementResultsUpdateType))

#include <mpi.h>
#include <functional>

namespace sgpp {
namespace datadriven {

enum MPI_COMMAND_ID {
  NULL_COMMAND,
  UPDATE_GRID,
  MERGE_GRID,
  ASSIGN_BATCH,
  COMPUTE_UPDATE_SYSTEM_MATRIX_DECOMPOSITION,
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
  size_t gridversion;
  size_t classIndex;
  size_t listLength;
  RefinementResultsUpdateType updateType;

  unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE];
};

struct RefinementResultCholeskyNetworkMessage {
  size_t matrixWidth;
  size_t matrixHeight;
  size_t offset;

  unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE - 3 * sizeof(size_t)];
};

struct MergeGridNetworkMessage {
  size_t gridversion;
  size_t classIndex;
  size_t payloadOffset;
  size_t payloadLength;
  size_t batchSize;
  size_t batchOffset;
  size_t alphaTotalSize;

  unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
      - 7 * sizeof(size_t))];
};

struct AssignBatchNetworkMessage {
  size_t batchOffset;
  size_t batchSize;
  bool doCrossValidation;
};

struct AssignSystemMatrixUpdateNetworkMessage {
  size_t classIndex;
  size_t gridversion;
};

static_assert(sizeof(MergeGridNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Merge Grid Network Message too long.");
static_assert(sizeof(RefinementResultNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Refinement result Network Message too long.");
static_assert(sizeof(RefinementResultCholeskyNetworkMessage) <= MPI_PACKET_MAX_PAYLOAD_SIZE,
              "Refinement result Cholesky Network Message too long.");
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_NETWORKMESSAGEDATA_H
