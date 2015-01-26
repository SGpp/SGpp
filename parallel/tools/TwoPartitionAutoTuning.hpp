/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TWOPARTIONAUTOTUNING_HPP
#define TWOPARTIONAUTOTUNING_HPP

#include <cstddef>

namespace sg {
  namespace parallel {
    class TwoPartitionAutoTuning {
      public:
        /**
         * Constructor
         *
         * @param problemSize contains the overall size which should be partitioned
         * @param partition2Divider the second partition divider, partition2's size is a multiple
         * @param retune_cycles number of iteration after which the problem's separation is re-considered
         */
        TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles);

        /**
         * Destructor
         */
        virtual ~TwoPartitionAutoTuning();

        /**
         * get problem size
         *
         * @return problem size to should partitioned
         */
        size_t getProblemSize();

        /**
         * sets the problem size
         *
         * @param problemSize problem size
         */
        void setProblemSize(size_t problemSize);

        /**
         * Update execution times in order to allow
         * a new calculation of the partition sizes
         *
         * @param timePartition1 time needed for partition 1
         * @param timePartition2 time needed for partition 2
         */
        void setExecutionTimes(double timePartition1, double timePartition2);

        /**
         * gets size of partition 1 based on the currently stored
         * runtimes for partition 1 and 2
         *
         * @return size of partition 1
         */
        size_t getPartition1Size();

        /**
         * resets all auto tuning parameter
         */
        virtual void resetAutoTuning();


        /**
         * set the possible divider of partition 2
         *
         * @param partition2Divider the divider of partition 2
         */
        void setPartition2Divider(size_t partition2Divider);
      protected:
        /**
         * @brief autoTune this function is called
         * - after the executionTimes are set to adapt the partitions or
         * - after the total size changed
         */
        virtual void autoTune() = 0;

        void doTune();
        /// size of partition1
        size_t _sizePartition1;

        /// store problemsize
        size_t _problemSize;

        /// time needed to execute partition 1
        double _timePartition1;

        /// time needed to execute partition 2
        double _timePartition2;

        /// store required divider of partition 2
        size_t _partition2Divider;

        /// counter for timer updates
        size_t _tuneCounter;

        /// number of updates that cause a tuning update
        size_t _retune;
    };

  }

}

#endif /* TWOPARTIONAUTOTUNING_HPP */
