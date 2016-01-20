// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DYNAMICTWOPARTITIONAUTOTUNING_HPP
#define DYNAMICTWOPARTITIONAUTOTUNING_HPP

#include <sgpp/parallel/tools/TwoPartitionAutoTuning.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    class DynamicTwoPartitionAutoTuning : public TwoPartitionAutoTuning {
      public:
        DynamicTwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles);
        virtual void resetAutoTuning();
      protected:
        void autoTune();
        double _partition2_speedup;
    };

  }
}
#endif // DYNAMICTWOPARTITIONAUTOTUNING_HPP