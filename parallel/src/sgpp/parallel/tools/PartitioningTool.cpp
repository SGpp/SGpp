// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/PartitioningTool.hpp>

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <sgpp/base/exception/operation_exception.hpp>
#include <iostream>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#ifdef USE_MPI
#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    PartitioningTool::PartitioningTool() {
    }

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
    void PartitioningTool::getPartitionSegment(size_t totalSize, size_t segmentCount, size_t segmentNumber, size_t* size, size_t* offset, size_t blockSize) {
      size_t end;
      getPartitionSegment(0, totalSize, segmentCount, segmentNumber, offset, &end, blockSize);
      *size = end - *offset;
    }

    void PartitioningTool::getPartitionSegment(size_t start, size_t end, size_t segmentCount, size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd, size_t blockSize) {
      size_t totalSize = end - start;

      // check for valid input
      if (blockSize == 0 ) {
        throw SGPP::base::operation_exception("blockSize must not be zero!");
      }

      if (totalSize % blockSize != 0 ) {
        //std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
        throw SGPP::base::operation_exception("totalSize must be divisible by blockSize without remainder, but it is not!");
      }

      // do all further calculations with complete blocks
      size_t blockCount = totalSize / blockSize;

      size_t blockSegmentSize = blockCount / segmentCount;
      size_t remainder = blockCount - blockSegmentSize * segmentCount;
      size_t blockSegmentOffset = 0;

      if (segmentNumber < remainder) {
        blockSegmentSize++;
        blockSegmentOffset = blockSegmentSize * segmentNumber;
      } else {
        blockSegmentOffset = remainder * (blockSegmentSize + 1) + (segmentNumber - remainder) * blockSegmentSize;
      }

      *segmentStart = start + blockSegmentOffset * blockSize;
      *segmentEnd = *segmentStart + blockSegmentSize * blockSize;
    }

    void PartitioningTool::getOpenMPPartitionSegment(size_t totalSize, size_t* size, size_t* offset, size_t blocksize) {
      size_t end;
      getOpenMPPartitionSegment(0, totalSize, offset, &end, blocksize);
      *size = end - *offset;
    }

    void PartitioningTool::getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
      size_t threadCount = 1;
      size_t myThreadNum = 0;
#ifdef _OPENMP
      threadCount = omp_get_num_threads();
      myThreadNum = omp_get_thread_num();
#endif
      getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
    }

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#ifdef USE_MPI
    void PartitioningTool::getMPIPartitionSegment(size_t totalSize, size_t* size, size_t* offset, size_t blocksize) {
      size_t end;
      getMPIPartitionSegment(0, totalSize, offset, &end, blocksize);
      *size = end - *offset;
    }

    void PartitioningTool::getMPIPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
      size_t myRank = 0;
      size_t numRanks = 1;

      myRank = SGPP::parallel::myGlobalMPIComm->getMyRank();
      numRanks = SGPP::parallel::myGlobalMPIComm->getNumRanks();

      getPartitionSegment(start, end, numRanks, myRank, segmentStart, segmentEnd, blocksize);
    }
#endif

    void PartitioningTool::calcDistribution(size_t totalSize, size_t numChunks, int* sizes, int* offsets, size_t blocksize) {
      for (size_t chunkID = 0; chunkID < numChunks; ++chunkID) {
        size_t size;
        size_t offset;
        getPartitionSegment(totalSize, numChunks, chunkID, &size, &offset, blocksize);
        sizes[chunkID] = (int)size;
        offsets[chunkID] = (int)offset;
      }
    }

#ifdef USE_MPI
    void PartitioningTool::calcMPIChunkedDistribution(size_t totalSize, size_t numChunksPerProc, int* sizes, int* offsets, size_t blocksize) {
      size_t numRanks = 1;

      numRanks = myGlobalMPIComm->getNumRanks();

      size_t procSize;
      size_t procOffset;

      for (size_t proc = 0; proc < numRanks; proc++) {
        getPartitionSegment(totalSize, numRanks, proc, &procSize, &procOffset, blocksize);
        calcDistribution(procSize, numChunksPerProc, &sizes[numChunksPerProc * proc], &offsets[numChunksPerProc * proc], blocksize);

        for (size_t i = 0; i < numChunksPerProc; i++) {
          offsets[numChunksPerProc * proc + i] += (int)(procOffset);
        }
      }
    }
#endif
  }
}