/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/tools/TwoPartitionAutoTuning.hpp>
#include <iostream>
#include <algorithm>
#include <iomanip>

#define INITIAL_SPEEDUP_PARTITION_2 20.0;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    TwoPartitionAutoTuning::TwoPartitionAutoTuning(size_t problemSize, size_t partition2Divider, size_t retune_cycles)
      : _problemSize(problemSize),
        _timePartition1(0.0),
        _timePartition2(0.0),
        _partition2Divider(partition2Divider),
        _tuneCounter(0),
        _retune(retune_cycles) {
    }

    TwoPartitionAutoTuning::~TwoPartitionAutoTuning() {
    }

    size_t TwoPartitionAutoTuning::getProblemSize() {
      return _problemSize;
    }

    void TwoPartitionAutoTuning::setProblemSize(size_t problemSize) {
      if (_problemSize == problemSize) {
        return;
      }

      _problemSize = problemSize;
      doTune();
    }

    void TwoPartitionAutoTuning::setPartition2Divider(size_t partition2Divider) {
      if (_partition2Divider == partition2Divider) {
        return;
      }

      _partition2Divider = partition2Divider;
      doTune();
    }

    void TwoPartitionAutoTuning::doTune() {
      autoTune();
      std::ios_base::fmtflags f = std::cout.flags();
      std::streamsize prec = std::cout.precision();
      std::cout << "AUTOTUNING-PARTITION-SIZES (" << std::setw(8) << _problemSize << "):"
                << std::setiosflags(std::ios::fixed) << std::setprecision(3)
                << " Time1: " << _timePartition1
                << " Size1: " << std::setw(8) << getPartition1Size()
                << std::setprecision(1)
                << "(" << 100.0 * (double)getPartition1Size() / (double)_problemSize << "%); ";
      std::cout << " Time2: " << std::setprecision(3) << _timePartition2
                << " Size2: " << std::setw(8) << _problemSize - getPartition1Size()
                << std::setprecision(1)
                << " (" << 100.0 * (double)(_problemSize - getPartition1Size()) / (double)_problemSize << "%)" << std::endl;
      std::cout.flags(f);
      std::cout.precision(prec);
    }

    void TwoPartitionAutoTuning::setExecutionTimes(double timePartition1, double timePartition2) {
      _timePartition1 += timePartition1;
      _timePartition2 += timePartition2;
      _tuneCounter++;

      if (_tuneCounter % _retune == 0) {
        doTune();
        _tuneCounter = 0;
        _timePartition1 = 0;
        _timePartition2 = 0;
      }
    }

    void TwoPartitionAutoTuning::resetAutoTuning() {
      _timePartition1 = 0.0;
      _timePartition2 = 0.0;
      _tuneCounter = 0;
      doTune();
    }

    size_t TwoPartitionAutoTuning::getPartition1Size() {
      return _sizePartition1;
    }

  }

}
