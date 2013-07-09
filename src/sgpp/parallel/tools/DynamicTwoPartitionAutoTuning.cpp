/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#define INITIAL_SPEEDUP_PARTITION_2 20.0

#include "DynamicTwoPartitionAutoTuning.hpp"
namespace sg {
  namespace parallel {

    DynamicTwoPartitionAutoTuning::DynamicTwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles) :
      TwoPartitionAutoTuning(problemSize, partition2Divider, retune_cycles),
      _partition2_speedup(INITIAL_SPEEDUP_PARTITION_2) {
      doTune();
    }

    void DynamicTwoPartitionAutoTuning::resetAutoTuning() {
      _partition2_speedup = INITIAL_SPEEDUP_PARTITION_2;
      TwoPartitionAutoTuning::resetAutoTuning();
    }

    void DynamicTwoPartitionAutoTuning::autoTune() {

      // don't use an accelerator for smallest problems
      if (4 * _partition2Divider > _problemSize) {
        _sizePartition1 = _problemSize;
        return;
      }

      size_t partition1 = 0;
      size_t partition2 = 0;

      if (((_tuneCounter % _retune) == 0)) {
        if (_timePartition1 != 0 && _timePartition2 != 0) {
          double partition1_element_time = static_cast<double>(_sizePartition1) / static_cast<double>(_timePartition1);
          double partition2_element_time = static_cast<double>(_problemSize - _sizePartition1) / static_cast<double>(_timePartition2);
          _partition2_speedup = partition2_element_time / partition1_element_time;
        }
      }

      double normalized_workingset = static_cast<double>(_problemSize) / (_partition2_speedup + 0.93);

      partition2 = static_cast<size_t>(normalized_workingset * _partition2_speedup);

      if (partition2 < _partition2Divider / 2) {   // don't use an accelerator for too small chunks
        _sizePartition1 = _problemSize;
        return;
      }

      if (partition2 == _problemSize) {
        partition2 -= _partition2Divider;
      }

      size_t partition2_remainder = partition2 % _partition2Divider;

      if (partition2 + (_partition2Divider - partition2_remainder) >= _problemSize) {
        partition2 -=  partition2_remainder;
      } else {
        partition2 +=  (_partition2Divider - partition2_remainder);
      }

      partition1 = _problemSize - partition2;

      _sizePartition1 = partition1;
    }

  }
}
