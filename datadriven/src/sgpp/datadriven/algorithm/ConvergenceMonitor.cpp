// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>

namespace sgpp {
namespace datadriven {

ConvergenceMonitor::ConvergenceMonitor(double pDeclineThreshold,
                                       size_t pBufferSize,
                                       size_t pMinRefInterval)
    : nextRefCnt(0),
      minRefInterval(pMinRefInterval),
      validErrorSum1(0.0),
      validErrorSum2(0.0),
      trainErrorSum1(0.0),
      trainErrorSum2(0.0),
      validDiff(1.0),
      trainDiff(1.0),
      declineThreshold(pDeclineThreshold),
      bufferSize(pBufferSize) {}

ConvergenceMonitor::~ConvergenceMonitor() {}

void ConvergenceMonitor::pushToBuffer(double currentValidError,
                                      double currentTrainError) {
  if (validErrorDeclineBuffer.size() >= bufferSize) {
    validErrorDeclineBuffer.pop_back();
    trainErrorDeclineBuffer.pop_back();
  }
  validErrorDeclineBuffer.push_front(currentValidError);
  trainErrorDeclineBuffer.push_front(currentTrainError);
}

bool ConvergenceMonitor::checkConvergence() {
  bool result = false;

  if (validErrorDeclineBuffer.size() == bufferSize) {
    validErrorSum1 = 0.0;
    validErrorSum2 = 0.0;
    trainErrorSum1 = 0.0;
    trainErrorSum2 = 0.0;
    size_t maxIdx = static_cast<size_t>(validErrorDeclineBuffer.size() / 2);
    for (size_t idx = 0; idx < maxIdx; idx++) {
      validErrorSum1 += validErrorDeclineBuffer[idx];
      validErrorSum2 += validErrorDeclineBuffer[idx + maxIdx];
      trainErrorSum1 += trainErrorDeclineBuffer[idx];
      trainErrorSum2 += trainErrorDeclineBuffer[idx + maxIdx];
    }

    validDiff = validErrorSum2 / static_cast<double>(maxIdx) -
                validErrorSum1 / static_cast<double>(maxIdx);
    trainDiff = trainErrorSum2 / static_cast<double>(maxIdx) -
                trainErrorSum1 / static_cast<double>(maxIdx);
  }

  if (validDiff >= 0.0) {
    if (validDiff < declineThreshold) {
      result = true;
    }
  } else if (validDiff < 0.0) {
    if (trainDiff < 0.0) {
      result = true;
    }
  }
  return result;
}

}  // namespace datadriven
}  // namespace sgpp

