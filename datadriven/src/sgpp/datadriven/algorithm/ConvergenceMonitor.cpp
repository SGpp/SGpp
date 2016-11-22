// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>

#include <iostream>
//#include <cmath>

ConvergenceMonitor::ConvergenceMonitor(double pDeclineThreshold,
                                       size_t pBufferSize, size_t pMinRefInterval) 
    : nextRefCnt(0),
      minRefInterval(pMinRefInterval),
      validErrorSum1(0.0),
      validErrorSum2(0.0),
      trainErrorSum1(0.0),
      trainErrorSum2(0.0),
      validRatio(1.0),
      trainRatio(1.0),
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
    size_t maxIdx = (size_t)(validErrorDeclineBuffer.size()/2);
    for (size_t idx = 0; idx < maxIdx; idx++) {
      validErrorSum1 += validErrorDeclineBuffer[idx];
      validErrorSum2 += validErrorDeclineBuffer[idx+maxIdx];
      trainErrorSum1 += trainErrorDeclineBuffer[idx];
      trainErrorSum2 += trainErrorDeclineBuffer[idx+maxIdx];
    }  

    validRatio = validErrorSum2/double(maxIdx) - validErrorSum1/double(maxIdx);
    trainRatio = trainErrorSum2/double(maxIdx) - trainErrorSum1/double(maxIdx);

    /*std::cout << "validErrorSum1: " << validErrorSum1 << std::endl;
    std::cout << "validErrorSum2: " << validErrorSum2 << std::endl;
    std::cout << "validRatio: " << validRatio << std::endl;
    std::cout << "trainRatio: " << trainRatio << std::endl;*/
  }

  if (validRatio >= 0.0) { 
    if (validRatio < declineThreshold) {
      //std::cout << "crit. 1 refinement" << std::endl;
      result = true;
    }
  }
  /*else if (validRatio < 0.0) {
    if (trainRatio < 0.0) {
      //std::cout << "crit. 2 refinement" << std::endl;
      result = true;
    }
  }*/
  return result;

}


