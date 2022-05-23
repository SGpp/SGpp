// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <omp.h>
#include <algorithm>
// #include <chrono>
// #include <thread>
#include <iostream>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalMPI/OperationMultiEvalMPI.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>

namespace sgpp {
namespace datadriven {

OperationMultiEvalMPI::OperationMultiEvalMPI(base::Grid& grid, base::DataMatrix& dataset,
                                             OperationMultipleEvalType nodeImplType,
                                             OperationMultipleEvalSubType nodeImplSubType,
                                             bool verbose)
    : OperationMultipleEval(grid, dataset),
      nodeImplType(nodeImplType),
      nodeImplSubType(nodeImplSubType),
      dim(grid.getDimension()),
      verbose(verbose),
      duration(-1.0) {
  // create the kernel specific data structures for the current grid
  this->prepare();
}

// OperationMultiEvalMPI::OperationMultiEvalMPI(base::DataMatrix& fakeGridLevel,
//                                             base::DataMatrix& fakeGridIndex, size_t dim,
//                                             base::DataMatrix& dataset,
//                                             OperationMultiEvalMPIType type, bool verbose)
//    : dataset(dataset),
//      level(fakeGridLevel),
//      index(fakeGridIndex),
//      type(type),
//      dim(dim),
//      verbose(verbose),
//      duration(-1.0) {
//  // create the kernel specific data structures for the current grid
//  this->prepare();
//}

OperationMultiEvalMPI::~OperationMultiEvalMPI() {}

void OperationMultiEvalMPI::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  double start = MPI_Wtime();

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<size_t> metaInformation(3);
  std::vector<size_t> workChunk(2);

  if (rank == 0) {
    size_t chunkSize = 5000;
    size_t currentChunkStart = 0;
    int ranksFinished = 0;

    // while there is still work to do or results to receive
    while (ranksFinished < size - 1) {
// receive result meta information from last iteration
// int messageAvailable = false;
// MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &messageAvailable, 0);
// if (!messageAvailable) {
//   //        std::cout << "waiting.. " << std::endl;
//   std::this_thread::sleep_for(std::chrono::milliseconds(5));
//   continue;
// }

#ifdef __NEC__
      MPI_Recv(&metaInformation[0], 3, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
#else
      MPI_Recv(metaInformation.data(), 3, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
#endif
      int sender = static_cast<int>(metaInformation[0]);
      size_t resultChunkStart = metaInformation[1];
      size_t resultChunkEnd = metaInformation[2];
      size_t resultChunkRange = resultChunkEnd - resultChunkStart;

      if (verbose) {
        std::cout << "rank = " << rank << ", sender: " << sender
                  << " resultStart: " << resultChunkStart << " resultEnd: " << resultChunkEnd
                  << std::endl;
      }

      // am I receiving any results?
      if (resultChunkRange > 0) {
        if (verbose) {
          std::cout << "rank = " << rank << ", sender: " << sender << ": reading result values..."
                    << std::endl;
        }
        std::vector<double> resultTemp(resultChunkRange);
#ifdef __NEC__
        MPI_Recv(&resultTemp[0], static_cast<int>(resultChunkRange), MPI_DOUBLE, sender, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
        MPI_Recv(resultTemp.data(), static_cast<int>(resultChunkRange), MPI_DOUBLE, sender, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
        // accumulate results from the other rank
        for (size_t i = 0; i < resultChunkRange; i++) {
          result[resultChunkStart + i] = resultTemp[i];
        }

        if (verbose) {
          std::cout << "rank = " << rank << ", sender: " << sender << ": done reading result values"
                    << std::endl;
        }
      } else {
        if (verbose) {
          std::cout << "rank = " << rank << ", sender: " << sender << ": no results to read!"
                    << std::endl;
        }
      }

      // is there any new work still available?
      if (currentChunkStart < dataset.getNrows()) {
        size_t currentChunkEnd = std::min(currentChunkStart + chunkSize, dataset.getNrows());
        workChunk[0] = currentChunkStart;
        workChunk[1] = currentChunkEnd;
        if (verbose) {
          std::cout << "rank = " << rank << ", sender: " << sender
                    << ": work package start: " << currentChunkStart << ", end: " << currentChunkEnd
                    << std::endl;
        }
#ifdef __NEC__
        MPI_Send(&workChunk[0], 2, MPI_UNSIGNED_LONG_LONG, sender, 0, MPI_COMM_WORLD);
#else
        MPI_Send(workChunk.data(), 2, MPI_UNSIGNED_LONG_LONG, sender, 0, MPI_COMM_WORLD);
#endif

        currentChunkStart = currentChunkEnd;
      } else {
        // send empty range to signal that the worker can finish
        workChunk[0] = currentChunkStart;
        workChunk[1] = currentChunkStart;
        if (verbose) {
          std::cout << "rank = " << rank << ", sender: " << sender
                    << ": finish pseudo work package start/end: " << currentChunkStart << std::endl;
        }

#ifdef __NEC__
        MPI_Send(&workChunk[0], 2, MPI_UNSIGNED_LONG_LONG, sender, 0, MPI_COMM_WORLD);
#else
        MPI_Send(workChunk.data(), 2, MPI_UNSIGNED_LONG_LONG, sender, 0, MPI_COMM_WORLD);
#endif

        ranksFinished += 1;
      }
    }
  }
  this->duration = MPI_Wtime() - start;
}

void OperationMultiEvalMPI::multSlave(sgpp::base::DataVector& alpha) {
  double start = MPI_Wtime();

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<size_t> metaInformation(3);
  std::vector<size_t> workChunk(2);
  sgpp::base::DataVector result;

  size_t chunkStart = 0;
  size_t chunkEnd = 0;
  size_t chunkRange = 0;

  while (true) {
    // transmit result from last iteration // size, chunkStart, chunkEnd, data[]
    metaInformation[0] = rank;
    metaInformation[1] = chunkStart;
    metaInformation[2] = chunkEnd;

    if (verbose) {
      std::cout << "rank = " << rank << " to master, sender: " << rank
                << ", resultStart: " << chunkStart << " resultEnd: " << chunkEnd << std::endl;
    }

// send meta information to master
#ifdef __NEC__
    MPI_Send(&metaInformation[0], 3, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
#else
    MPI_Send(metaInformation.data(), 3, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
#endif
    // send the last result except for the initialization
    if (chunkRange > 0) {
      if (verbose) {
        std::cout << "rank = " << rank << " to master, sending result data..." << std::endl;
      }

      MPI_Send(result.getPointer(), static_cast<int>(chunkRange), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      if (verbose) {
        std::cout << "rank = " << rank << " to master, done sending result data!" << std::endl;
      }
    } else {
      if (verbose) {
        std::cout << "rank = " << rank << " to master, no result data to send!" << std::endl;
      }
    }

#ifdef __NEC__
    MPI_Recv(&workChunk[0], 2, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
    MPI_Recv(workChunk.data(), 2, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
    chunkStart = workChunk[0];
    chunkEnd = workChunk[1];
    chunkRange = chunkEnd - chunkStart;

    if (verbose) {
      std::cout << "rank = " << rank << " from master, work chunkStart: " << chunkStart
                << ", chunkEnd: " << chunkEnd << std::endl;
    }

    // empty chunk received, finishing...
    if (chunkRange == 0) {
      if (verbose) {
        std::cout << "rank = " << rank << ", shutting down" << std::endl;
      }
      break;
    }

    if (verbose) {
      std::cout << "rank = " << rank << ", doing calculations..." << std::endl;
    }

    // filter the dataset according to the range to work on
    base::DataMatrix datasetChunk(chunkRange, dataset.getNcols());
    for (size_t i = 0; i < chunkRange; i++) {
      for (size_t d = 0; d < dim; d++) {
        datasetChunk.set(i, d, dataset.get(chunkStart + i, d));
      }
    }

    // resize the result vector to matchthe chunk
    result.resize(chunkRange);

    // create appropriate node level multi eval implementation
    sgpp::base::OperationMultipleEval* nodeMultiEval;
    if (nodeImplType == OperationMultipleEvalType::STREAMING &&
        nodeImplSubType == OperationMultipleEvalSubType::DEFAULT) {
      //      nodeMultiEval = new datadriven::OperationMultiEvalStreaming(level, index, dim,
      //      datasetChunk);
      nodeMultiEval = new datadriven::OperationMultiEvalStreaming(grid, datasetChunk);
    } else {
      throw base::not_implemented_exception();
    }

    // calculate the result for the received range
    nodeMultiEval->mult(alpha, result);

    delete nodeMultiEval;

    if (verbose) {
      std::cout << "rank = " << rank << ", calculations finished" << std::endl;
    }
  }

  if (verbose) {
    std::cout << "rank = " << rank << ", master finished" << std::endl;
  }

  this->duration = MPI_Wtime() - start;
}

void OperationMultiEvalMPI::multTranspose(sgpp::base::DataVector& source,
                                          sgpp::base::DataVector& result) {
  double start = MPI_Wtime();

  //  base::DataMatrix levelChunk(level);
  //  base::DataMatrix indexChunk(level);
  //
  //  // create appropriate node level multi eval implementation
  //  sgpp::base::OperationMultiEval* nodeMultiEval;
  //  if (type == OperationMultiEvalMPIType::STREAMING) {
  //    nodeMultiEval =
  //        new datadriven::OperationMultiEvalStreaming(levelChunk, indexChunk, dim, dataset);
  //  }
  //  nodeMultiEval->multTranspose(source, result);
  //
  //  delete nodeMultiEval;

  throw base::not_implemented_exception();
  this->duration = MPI_Wtime() - start;
}

double OperationMultiEvalMPI::getDuration() { return this->duration; }

void OperationMultiEvalMPI::prepare() {}
}  // namespace datadriven
}  // namespace sgpp
