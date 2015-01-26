/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#include "StaticTwoPartitionAutoTuning.hpp"
#include <algorithm>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    StaticTwoPartitionAutoTuning::StaticTwoPartitionAutoTuning(size_t problemSize, double percentPartion1, size_t partition2Divider, size_t OutputFreq):
      TwoPartitionAutoTuning(problemSize, partition2Divider, OutputFreq),
      _percentPartion1(percentPartion1) {
      doTune();
    }

    void StaticTwoPartitionAutoTuning::autoTune() {
      size_t partition1 = (size_t)std::min<double>(((double)_problemSize) * _percentPartion1, (double)_problemSize);
      size_t partition2 = _problemSize - partition1;

      size_t partition2_remainder = partition2 % _partition2Divider;
      partition2 -=  partition2_remainder;
      partition1 = _problemSize - partition2;

      _sizePartition1 = partition1;
    }



  }
}
