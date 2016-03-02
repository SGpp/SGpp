// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/MPICommunicator.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <string>
namespace SGPP {
namespace parallel {

MPICommunicator::MPICommunicator(int myid, int ranks) : myid_(myid), ranks_(ranks) {}

MPICommunicator::~MPICommunicator() {}

void MPICommunicator::broadcastGridCoefficientsFromRank0(SGPP::base::DataVector& alpha) {
  MPI_Bcast(reinterpret_cast<void*>(alpha.getPointer()), static_cast<int>(alpha.getSize()),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPICommunicator::broadcastSPGridCoefficientsFromRank0(base::DataVectorSP& alpha) {
  MPI_Bcast(reinterpret_cast<void*>(alpha.getPointer()), static_cast<int>(alpha.getSize()),
            MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void MPICommunicator::reduceGridCoefficientsOnRank0(SGPP::base::DataVector& alpha) {
  if (myid_ == 0) {
    MPI_Reduce(MPI_IN_PLACE, reinterpret_cast<void*>(alpha.getPointer()),
               static_cast<int>(alpha.getSize()), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(reinterpret_cast<void*>(alpha.getPointer()), NULL, static_cast<int>(alpha.getSize()),
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

void MPICommunicator::reduceGridCoefficients(SGPP::base::DataVector& alpha) {
  MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<void*>(alpha.getPointer()),
                static_cast<int>(alpha.getSize()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MPICommunicator::sendGrid(std::string& serialized_grid, int dest_rank) {
  int nChars = static_cast<int>(serialized_grid.length());
  MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
  MPI_Send(reinterpret_cast<void*>(const_cast<char*>(serialized_grid.c_str())),
           static_cast<int>(serialized_grid.length()), MPI_CHAR, dest_rank, this->myid_,
           MPI_COMM_WORLD);
}

void MPICommunicator::broadcastGrid(std::string& serialized_grid) {
  int nChars = static_cast<int>(serialized_grid.length());

  for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++) {
    MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
    MPI_Send(reinterpret_cast<void*>(const_cast<char*>(serialized_grid.c_str())),
             static_cast<int>(serialized_grid.length()), MPI_CHAR, dest_rank, this->myid_,
             MPI_COMM_WORLD);
  }
}

void MPICommunicator::receiveGrid(std::string& serialized_grid) {
  char* recv_buffer;
  int count;
  MPI_Status status;

  MPI_Recv(reinterpret_cast<void*>((&count)), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
           MPI_COMM_WORLD, &status);
  recv_buffer = new char[count];
  MPI_Recv(reinterpret_cast<void*>(recv_buffer), count, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
           MPI_COMM_WORLD, &status);

  serialized_grid = "";
  serialized_grid.assign(recv_buffer);

  delete recv_buffer;
}

void MPICommunicator::broadcastGridStorage(std::string& serialized_grid_storage) {
  int nChars = static_cast<int>(serialized_grid_storage.length());

  for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++) {
    MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
    MPI_Send(reinterpret_cast<void*>(const_cast<char*>(serialized_grid_storage.c_str())),
             static_cast<int>(serialized_grid_storage.length()), MPI_CHAR, dest_rank, this->myid_,
             MPI_COMM_WORLD);
  }
}

void MPICommunicator::receiveGridStorage(std::string& serialized_grid_storage) {
  char* recv_buffer;
  int count;
  MPI_Status status;

  MPI_Recv(reinterpret_cast<void*>((&count)), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
           MPI_COMM_WORLD, &status);
  recv_buffer = new char[count];
  MPI_Recv(reinterpret_cast<void*>(recv_buffer), count, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
           MPI_COMM_WORLD, &status);

  serialized_grid_storage = "";
  serialized_grid_storage.assign(recv_buffer);

  delete recv_buffer;
}

void MPICommunicator::Barrier() { MPI_Barrier(MPI_COMM_WORLD); }

void MPICommunicator::Abort() { MPI_Abort(MPI_COMM_WORLD, -1); }

void MPICommunicator::broadcastControlFromRank0(char* ctrl) {
  MPI_Bcast(reinterpret_cast<void*>(ctrl), 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void MPICommunicator::dataVectorAllToAll(base::DataVector& alpha, int* distributionOffsets,
                                         int* distributionSizes) {
  size_t myRank = getMyRank();
  int mySendSize = distributionSizes[myRank];
  int mySendOffset = distributionOffsets[myRank];

  SGPP::base::DataVector tmp(alpha.getSize());
  double* sendbuf = alpha.getPointer();
  MPI_Allgatherv(&(sendbuf[mySendOffset]), mySendSize, MPI_DOUBLE, tmp.getPointer(),
                 distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);
  alpha.copyFrom(tmp);
}

void MPICommunicator::dataVectorSPAllToAll(base::DataVectorSP& alpha, int* distributionOffsets,
                                           int* distributionSizes) {
  size_t myRank = getMyRank();
  int mySendSize = distributionSizes[myRank];
  int mySendOffset = distributionOffsets[myRank];

  SGPP::base::DataVectorSP tmp(alpha.getSize());
  float* sendbuf = alpha.getPointer();
  MPI_Allgatherv(&(sendbuf[mySendOffset]), mySendSize, MPI_FLOAT, tmp.getPointer(),
                 distributionSizes, distributionOffsets, MPI_FLOAT, MPI_COMM_WORLD);
  alpha.copyFrom(tmp);
}

void MPICommunicator::IsendToAll(double* ptr, size_t size, int tag, MPI_Request* reqs) {
  for (size_t rank = 0; rank < getNumRanks(); rank++) {
    if (rank == getMyRank()) {
      reqs[rank] = MPI_REQUEST_NULL;
      continue;
    }

    int sizeAsInt = static_cast<int>(size);
    MPI_Isend(ptr, sizeAsInt, MPI_DOUBLE, static_cast<int>(rank), tag, MPI_COMM_WORLD, &reqs[rank]);
  }
}

void MPICommunicator::IrecvFromAll(double* ptr, size_t chunkSizePerProc, int* sizes, int* offsets,
                                   int* tag, MPI_Request* reqs) {
  for (size_t rank = 0; rank < getNumRanks(); rank++) {
    for (size_t i = 0; i < chunkSizePerProc; i++) {
      size_t reqIdx = rank * chunkSizePerProc + i;

      if (rank == getMyRank()) {
        reqs[reqIdx] = MPI_REQUEST_NULL;
      } else {
        MPI_Irecv(&ptr[offsets[reqIdx]], sizes[reqIdx], MPI_DOUBLE, static_cast<int>(rank),
                  tag[reqIdx], MPI_COMM_WORLD, &reqs[reqIdx]);
      }
    }
  }
}

void MPICommunicator::IsendToAllSP(float* ptr, size_t size, int tag, MPI_Request* reqs) {
  for (size_t rank = 0; rank < getNumRanks(); rank++) {
    if (rank == getMyRank()) {
      reqs[rank] = MPI_REQUEST_NULL;
      continue;
    }

    int sizeAsInt = static_cast<int>(size);
    MPI_Isend(ptr, sizeAsInt, MPI_FLOAT, static_cast<int>(rank), tag, MPI_COMM_WORLD, &reqs[rank]);
  }
}

void MPICommunicator::IrecvFromAllSP(float* ptr, size_t chunkSizePerProc, int* sizes, int* offsets,
                                     int* tag, MPI_Request* reqs) {
  for (size_t rank = 0; rank < getNumRanks(); rank++) {
    for (size_t i = 0; i < chunkSizePerProc; i++) {
      size_t reqIdx = rank * chunkSizePerProc + i;

      if (rank == getMyRank()) {
        reqs[reqIdx] = MPI_REQUEST_NULL;
      } else {
        MPI_Irecv(&ptr[offsets[reqIdx]], sizes[reqIdx], MPI_FLOAT, static_cast<int>(rank),
                  tag[reqIdx], MPI_COMM_WORLD, &reqs[reqIdx]);
      }
    }
  }
}

void MPICommunicator::putToAll(double* ptr, size_t winOffset, size_t count, MPI_Win win) {
  for (size_t i = 0; i < getNumRanks(); i++) {
    int countAsInt = static_cast<int>(count);
    int winOffsetAsInt = static_cast<int>(winOffset);
    MPI_Put(ptr, countAsInt, MPI_DOUBLE, static_cast<int>(i), winOffsetAsInt, countAsInt,
            MPI_DOUBLE, win);
  }
}

void MPICommunicator::putToAllSP(float* ptr, size_t winOffset, size_t count, MPI_Win win) {
  for (size_t i = 0; i < getNumRanks(); i++) {
    int countAsInt = static_cast<int>(count);
    int winOffsetAsInt = static_cast<int>(winOffset);
    MPI_Put(ptr, countAsInt, MPI_FLOAT, static_cast<int>(i), winOffsetAsInt, countAsInt, MPI_FLOAT,
            win);
  }
}

size_t MPICommunicator::getMyRank() { return this->myid_; }

size_t MPICommunicator::getNumRanks() { return this->ranks_; }

void MPICommunicator::waitForAllRequests(size_t size, MPI_Request* reqs) {
  if (MPI_Waitall(static_cast<int>(size), reqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
    std::cout << "communication error in waitall" << std::endl;
    throw SGPP::base::operation_exception("Communication Error");
  }
}

void MPICommunicator::waitForAnyRequest(size_t size, MPI_Request* reqs, int* result) {
  MPI_Waitany(static_cast<int>(size), reqs, result, MPI_STATUS_IGNORE);
}

void MPICommunicator::allreduceSum(base::DataVector& source, base::DataVector& result) {
  if (source.getSize() != result.getSize()) {
    std::cout << "DataVector sizes do not match in allreduce!" << std::endl;
    throw SGPP::base::operation_exception("DataVector sizes do not match in allreduce!");
  }

  MPI_Allreduce(source.getPointer(), result.getPointer(), static_cast<int>(source.getSize()),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MPICommunicator::allreduceSumSP(base::DataVectorSP& source, base::DataVectorSP& result) {
  if (source.getSize() != result.getSize()) {
    std::cout << "DataVector sizes do not match in allreduce!" << std::endl;
    throw SGPP::base::operation_exception("DataVector sizes do not match in allreduce!");
  }

  MPI_Allreduce(source.getPointer(), result.getPointer(), static_cast<int>(source.getSize()),
                MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}
}  // namespace parallel
}  // namespace SGPP
