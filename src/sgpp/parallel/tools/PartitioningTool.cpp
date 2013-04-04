/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifdef USE_MPI
#include "parallel/tools/MPI/SGppMPITools.hpp"
#endif
#include "parallel/tools/PartitioningTool.hpp"
#include "base/exception/operation_exception.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <iostream>

namespace sg {
  namespace parallel {

    PartitioningTool::PartitioningTool() {
    }

    void PartitioningTool::getPartitionSegment(int totalSize, int segmentCount, int segmentNumber, size_t* size, size_t* offset, size_t blockSize) {
      size_t end;
      getPartitionSegment(0, totalSize, segmentCount, segmentNumber, offset, &end, blockSize);
      *size = end - *offset;
    }

    void PartitioningTool::getPartitionSegment(int start, int end, int segmentCount, int segmentNumber, size_t* segmentStart, size_t* segmentEnd, size_t blockSize) {
      int totalSize = end - start;

      // check for valid input
      if (totalSize % blockSize != 0 ) {
        std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
        throw sg::base::operation_exception("totalSize must be divisible by blockSize without remainder, but it is not!");
      }

      if (blockSize == 0 ) {
        throw sg::base::operation_exception("blockSize must not be zero!");
      }

      // do all further calculations with complete blocks
      int blockCount = totalSize / (int)blockSize;

      int blockSegmentSize = blockCount / segmentCount;
      int remainder = blockCount - blockSegmentSize * segmentCount;
      int blockSegmentOffset = 0;

      if (segmentNumber < remainder) {
        blockSegmentSize++;
        blockSegmentOffset = blockSegmentSize * segmentNumber;
      } else {
        blockSegmentOffset = remainder * (blockSegmentSize + 1) + (segmentNumber - remainder) * blockSegmentSize;
      }

      *segmentStart = start + blockSegmentOffset * blockSize;
      *segmentEnd = *segmentStart + blockSegmentSize * blockSize;
    }

    void PartitioningTool::getOpenMPPartitionSegment(int totalSize, size_t* size, size_t* offset, size_t blocksize) {
      size_t end;
      getOpenMPPartitionSegment(0, totalSize, offset, &end, blocksize);
      *size = end - *offset;
    }

    void PartitioningTool::getOpenMPPartitionSegment(int start, int end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
      int threadCount = 1;
      int myThreadNum = 0;
#ifdef _OPENMP
      threadCount = omp_get_num_threads();
      myThreadNum = omp_get_thread_num();
#endif
      getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart, segmentEnd, blocksize);
    }

#ifdef USE_MPI
    void PartitioningTool::getMPIPartitionSegment(int totalSize, size_t* size, size_t* offset, size_t blocksize) {
      size_t end;
      getMPIPartitionSegment(0, totalSize, offset, &end, blocksize);
      *size = end - *offset;
    }

    void PartitioningTool::getMPIPartitionSegment(int start, int end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
      int myRank = 0;
      int numRanks = 1;
#ifdef USE_MPI
      myRank = sg::parallel::myGlobalMPIComm->getMyRank();
      numRanks = sg::parallel::myGlobalMPIComm->getNumRanks();
#endif
      getPartitionSegment(start, end, numRanks, myRank, segmentStart, segmentEnd, blocksize);
    }

    void PartitioningTool::calcDistribution(int totalSize, int numChunks, int* sizes, int* offsets, size_t blocksize) {
      for (int chunkID = 0; chunkID < numChunks; ++chunkID) {
        size_t size;
        size_t offset;
        getPartitionSegment(totalSize, numChunks, chunkID, &size, &offset, blocksize);
        sizes[chunkID] = (int)size;
        offsets[chunkID] = (int)offset;
      }
    }

    void PartitioningTool::calcMPIChunkedDistribution(int totalSize, int numChunksPerProc, int* sizes, int* offsets, size_t blocksize) {
      int numRanks = 1;
#ifdef USE_MPI
      numRanks = myGlobalMPIComm->getNumRanks();
#endif
      int sizesGlobal[numRanks];
      int offsetsGlobal[numRanks];
      calcMPIChunkedDistribution(totalSize, numChunksPerProc, sizes, offsets, sizesGlobal, offsetsGlobal, blocksize);
    }

    void PartitioningTool::calcMPIChunkedDistribution(int totalSize, int numChunksPerProc, int* sizes, int* offsets, int* sizesGlobal, int* offsetsGlobal, size_t blocksize) {
      int numRanks = 1;
#ifdef USE_MPI
      numRanks = myGlobalMPIComm->getNumRanks();
#endif
      calcDistribution(totalSize, numRanks, sizesGlobal, offsetsGlobal, blocksize);
      int currentGlobalOffset = 0;

      for (int proc = 0; proc < numRanks; proc++) {
        calcDistribution(sizesGlobal[proc], numChunksPerProc, &sizes[numChunksPerProc * proc], &offsets[numChunksPerProc * proc], blocksize);
        int thisProcSize = 0;

        for (int i = 0; i < numChunksPerProc; i++) {
          offsets[numChunksPerProc * proc + i] += currentGlobalOffset;
          thisProcSize += sizes[numChunksPerProc * proc + i];
        }

        sizesGlobal[proc] = numChunksPerProc;
        offsetsGlobal[proc] = proc * numChunksPerProc;
        currentGlobalOffset += thisProcSize;
      }
    }
#endif

  }
}
