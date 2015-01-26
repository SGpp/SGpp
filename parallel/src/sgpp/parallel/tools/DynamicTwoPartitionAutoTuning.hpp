/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DYNAMICTWOPARTITIONAUTOTUNING_HPP
#define DYNAMICTWOPARTITIONAUTOTUNING_HPP

#include "TwoPartitionAutoTuning.hpp"

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
