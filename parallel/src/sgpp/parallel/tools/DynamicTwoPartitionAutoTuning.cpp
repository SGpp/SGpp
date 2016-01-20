// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define INITIAL_SPEEDUP_PARTITION_2 20.0

#include <sgpp/parallel/tools/DynamicTwoPartitionAutoTuning.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {
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

        // If for some reason partition2 doesn't get any chunks (too small problemsize e.g.),
        // _partition2_speedup will be 0. If something (problem size, workload per element, etc.)
        // changes, speedup will stay at 0 because the size of partition 2 will
        // stay at 0 and therefore partition2_element_time == 0.
        // Thus, if the time for partition1 is 3 times bigger than the one for partition2, we
        // can try again with a small chunksize for partition2. The speedup assigned inside the if results in a size of
        // artition 2 that is near 2*_partition2Divider
        // Then, autotuning can continue normally, as there are valid partition2_element_time values.
        if (_partition2_speedup == 0 && _timePartition1 / _timePartition2 > 3) {
          _partition2_speedup = 2.0 * static_cast<double>(_partition2Divider) / static_cast<double>(_problemSize - _partition2Divider);
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