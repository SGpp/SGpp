// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STATICTWOPARTITIONAUTOTUNING_HPP
#define STATICTWOPARTITIONAUTOTUNING_HPP

#include <sgpp/parallel/tools/TwoPartitionAutoTuning.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {


    class StaticTwoPartitionAutoTuning : public TwoPartitionAutoTuning {

        /**
         * Constructor for static load balancing
         *
         * @param problemSize contains the overall size which should be partitioned
         * @param percentPartion1 how big is the first, non accelerated portion, values must be from 0.0 to 1.0?
         * @param partition2Divider the second partition divider, partition2's size is a multiple
         * @param OutputFreq how often should we print timings?
         */
        StaticTwoPartitionAutoTuning(size_t problemSize, double percentPartion1, size_t partition2Divider, size_t OutputFreq);

      protected:
        void autoTune();
        /// static, percent threshold of partition 1 (values: 0.0 to 1.0)
        double _percentPartion1;
    };
  }
}

#endif // STATICTWOPARTITIONAUTOTUNING_HPP