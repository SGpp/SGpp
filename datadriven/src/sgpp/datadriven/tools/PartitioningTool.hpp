/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef PARTITIONINGTOOL_HPP
#define PARTITIONINGTOOL_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * The methods in this class calculate size and offset of a segment for a partition of a domain.
     *
     * The domain can be either specified by its size (@a totalSize) or by start (including) and end (excluding) indexes.
     * Then, the number of resulting segments (@a segmentCount) and the number of the desired segment (@a segmentNumber)
     * have to be passed. The results are stored into @a size and @a offset. The last (and optional) parameter blocksize
     * specifies how the segments should be aligned to blocks. The complete size must be evenly divisible by the blocksize.
     *
     * Segments are distributed as equally as possible (the difference between the minimum and maximum number of
     * items is at most @a blocksize). When just dividing (integer division) and leaving the rest to the last segment, this
     * segment could have twice as much to do as all the others (example: totalSize=127, segmentCount = 16, blockSize = 1)
     *
     * @param totalSize size of domain that's to be partitioned
     * @param start start of domain that's to be partitioned, including
     * @param end end of domain that's to be partitioned, excluding
     * @param segmentCount number of segments to partition the domain into
     * @param segmentNumber specifies the number of the fragment for which to calculate size and offset
     * @param size output variable to put resulting size into
     * @param offset output variable to put resulting offset into
     */
    class PartitioningTool {
      public:
        PartitioningTool();
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
        static void getPartitionSegment(size_t totalSize, size_t segmentCount, size_t segmentNumber, size_t* size, size_t* offset, size_t blocksize = 1);
        static void getPartitionSegment(size_t start, size_t end, size_t segmentCount, size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd, size_t blocksize = 1);

        /**
         * @brief getOpenMPLoopPartitionSegment uses the number of OpenMP Threads and the threads id for segmentCount and segmentNumber
         *
         * Call this function inside a parallel openmp region. Can also be called if only one Thread is active or even if OpenMP is disabled,
         * then the result is one single partition.
         *
         */
        static void getOpenMPPartitionSegment(size_t totalSize, size_t* size, size_t* offset, size_t blocksize = 1);
        static void getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize = 1);
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
        /**
         * @brief calcDistribution calculates a distribution of a domain of size @a totalSize into @a numCunks chunks and
         * fills the two arrays @a sizes and @a offsets with the respective offsets and sizes (both arrays have to be
         * already allocated and must be of size @a numChunks). When blocksize greater than 1, then the resulting sizes
         * are a multiple of this blocksize.
         *
         * Example:
         *
         * |. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|
         * calcDistribution(50, 3, sizes, offsets, 5) ->
         * |. . . . .-. . . . .-. . . . .-. . . . .|. . . . .-. . . . .-. . . . .|. . . . .-. . . . .-. . . . .|
         *
         * sizes: [20,15,15]
         * offsets: [0,20,35]
         *
         * @param totalSize size of domain to distribute
         * @param numChunks
         * @param sizes output array to store resulting distribution sizes (array size must be numChunks)
         * @param offsets output array to store resulting distribution offsets (array size must be numChunks)
         * @param blocksize resulting sizes are a multiple of this blocksize.
         */
        static void calcDistribution(size_t totalSize, size_t numChunks, int* sizes, int* offsets, size_t blocksize = 1);
#ifdef USE_MPI
        static void calcMPIChunkedDistribution(size_t totalSize, size_t numChunksPerProc, int* sizes, int* offsets, size_t blocksize);

        /**
         * @brief getMPIPartitionSegment uses the number of MPI processes and the MPI rank for segmentCount and segmentNumber
         *
         * This function can also be used if MPI is disabled, then the result is one single partition.
         */
        static void getMPIPartitionSegment(size_t totalSize, size_t* size, size_t* offset, size_t blocksize = 1);
        static void getMPIPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd, size_t blocksize = 1);

#endif

        static void calcAlmostBlockedDistribution(size_t totalSize, size_t numChunksPerProc, int* sizes, int* offsets, size_t blocksize);
    };

  }

}

#endif // PARTITIONINGTOOL_HPP
