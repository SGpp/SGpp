// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONVERGENCEMONITOR_HPP
#define CONVERGENCEMONITOR_HPP

#include <deque>
#include <sgpp/globaldef.hpp>

/**
 * A monitor to decide if a learning algorithm has converged. The
 * convergence criterion is based on the comparison of error
 * measurements throughout the training process.
 */
class ConvergenceMonitor {
  public:
   /**
    * Constructor.
    *
    * @param pDeclineThreshold The convergence threshold
    * @param pBufferSize Number of error measurements which are considered 
    *        for convergence check
    * @param pMinRefInterval Minimum number of iterations before next refinement 
    *        is allowed to be performed
    */
    ConvergenceMonitor(double pDeclineThreshold,
                       size_t pBufferSize, size_t pMinRefInterval);
   /**
    * Destructor.
    */
    virtual ~ConvergenceMonitor();
  
   /**
    * Stores the current error values in the buffer. If the buffer 
    * has reached the maximum size, the oldest values are removed.
    * 
    * @param
    * @param
    */
    void pushToBuffer(double currentValidError,
                      double currentTrainError);
   /**
    * Examines the convergence criterion with the current
    * error observations.
    *
    * @return True if converged, false otherwise
    */
    bool checkConvergence();

    // counts how many iterations are yet to be performed until
    // next refinement can be triggered (only required if minRefInterval > 0 is chosen)
    size_t nextRefCnt;
    size_t minRefInterval;
    // stores the latest validation error observations
    std::deque<double> validErrorDeclineBuffer;
    // stores the latest training error observations
    std::deque<double> trainErrorDeclineBuffer;

  private:
    // old validation error
    double validErrorSum1;
    // new validation error
    double validErrorSum2;
    // old training error
    double trainErrorSum1;
    // new training error
    double trainErrorSum2;
    // difference between validation error measurements
    double validDiff;
    // difference between training error measurements
    double trainDiff;

    double declineThreshold;
    size_t bufferSize;   
	
};

#endif /* CONVERGENCEMONITOR_HPP */
