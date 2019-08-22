// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#define MPI_PACKET_MAX_PAYLOAD_SIZE 4096
#define MPI_MASTER_RANK 0
#define MPI_MAX_PROCESSOR_NAME_LENGTH 256
#define MPI_TAG_HIGH_PRIORITY_NO_BLOCK 42
#define MPI_TAG_STANDARD_COMMAND 41

#define REFINENEMT_RESULT_PAYLOAD_SIZE (MPI_PACKET_MAX_PAYLOAD_SIZE\
                                   - 3 * sizeof(size_t)\
                                   - sizeof(RefinementResultsUpdateType))
#include <sgpp/globaldef.hpp>

#include <mpi.h>
#include <functional>

namespace sgpp {
namespace datadriven {

/**
 * Different commands sent over MPI to allow the receiver to identify the message's contents.
 */
enum MPI_COMMAND_ID {
  /**
   * Used to identify packets where the command id has not been set
   */
      nullptr_COMMAND,
  /**
   * A packet that contains changes to the grid or the system matrix
   */
      UPDATE_GRID,
  /**
   * A packet that contains results from a training.
   */
      MERGE_GRID,
  /**
   * A packet that assigns a batch with specified parameters to a worker
   */
      ASSIGN_BATCH,
  /**
   * A packet that assigns updating the system matrix decomposition to a worker
   */
      COMPUTE_UPDATE_SYSTEM_MATRIX_DECOMPOSITION,
  /**
   * A packet that instructs a worker to shutdown its MPI facilities after completing all requests.
   */
      SHUTDOWN,
  /**
   * A confirmation packet sent by a worker to acknowledge all requests were completed.
   */
      WORKER_SHUTDOWN_SUCCESS
};

/**
 * A packet sent over MPI, using a command as a descriptor, and a wrapped package in the payload
 * for data.
 */
struct MPI_Packet {
  /**
   * The MPI command of this specific packet.
   */
  MPI_COMMAND_ID commandID;
  /**
   * The packet's data segment.
   */
  unsigned char payload[MPI_PACKET_MAX_PAYLOAD_SIZE];
};

/**
 * The type of message received in a UPDATE_GRID message type.
 */
enum RefinementResultsUpdateType {
  /**
   * Packet contains a set of coordinates for newly created grid points.
   */
      ADDED_GRID_POINTS_LIST,
  /**
   * Packet contains a set of indices for grid points which were deleted.
   */
      DELETED_GRID_POINTS_LIST,
  /**
   * Packet contains a part of an updated system matrix decomposition.
   */
      SYSTEM_MATRIX_DECOMPOSITION
};

/**
 * Packet wrapped in an UPDATE_GRID MPI_Packet, containing segmented changes for a specified class.
 */
struct RefinementResultNetworkMessage {
  /**
   * The version of the grid to set after applying the results. If more segments follow,
   * this version will be set to inconsistent.
   */
  size_t gridversion;
  /**
   * The index of the class to update the grid for.
   */
  size_t classIndex;
  /**
   * The number of changes contained in this specific packet.
   */
  size_t listLength;
  /**
   * The type of changes contained in this specific packet.
   */
  RefinementResultsUpdateType updateType;

  /**
   * The packet's data segment.
   */
  unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE];
};

/**
 * Packet wrapped in a RefinementResultNetwork Message that contains additional
 * information required when updating the system matrix.
 */
struct RefinementResultSystemMatrixNetworkMessage {
  /**
   * The new target system matrix width.
   */
  size_t matrixWidth;
  /**
   * The new target system matrix height.
   */
  size_t matrixHeight;
  /**
   * The offset from the start of the matrix to copy to.
   */
  size_t offset;

  /**
   * The data to be copied.
   */
  unsigned char payload[REFINENEMT_RESULT_PAYLOAD_SIZE - 3 * sizeof(size_t)];
};

/**
 * Packet wrapper in MPI_Packet containing segmented data from the alpha vector of the trained
 * system.
 */
struct MergeGridNetworkMessage {
  /**
   * The version of the grid where the training took place.
   */
  size_t gridversion;
  /**
   * The index of the class that was trained.
   */
  size_t classIndex;
  /**
   * The offset from the start of the alpha vector where this segment of data starts.
   */
  size_t payloadOffset;
  /**
   * The number of values contained in this segment.
   */
  size_t payloadLength;
  /**
   * The size of the batch that was trained with.
   */
  size_t batchSize;
  /**
   * The offset of the batch that was trained with.
   */
  size_t batchOffset;
  /**
   * The total size of the alpha vector over all segments.
   */
  size_t alphaTotalSize;

  /**
   * The packet's data.
   */
  unsigned char payload[(MPI_PACKET_MAX_PAYLOAD_SIZE
      - 7 * sizeof(size_t))];
};

/**
 * Message wrapped in MPI_Packet specifying an order to a worker to train from a batch.
 */
struct AssignBatchNetworkMessage {
  /**
   * The offset from the start of the training set to learn from.
   */
  size_t batchOffset;
  /**
   * The size of the batch to learn from.
   */
  size_t batchSize;
  /**
   * Whether to do cross validation.
   */
  bool doCrossValidation;
};

/**
 * Message wrapped in MPI_Packet specifying an order to a worker to update a class' system matrix
 * decomposition
 */
struct AssignSystemMatrixUpdateNetworkMessage {
  /**
   * The class to execute the system matrix decomposition update for.
   */
  size_t classIndex;
  /**
   * The new grid version to set after the update.
   */
  size_t gridversion;
};

}  // namespace datadriven
}  // namespace sgpp
