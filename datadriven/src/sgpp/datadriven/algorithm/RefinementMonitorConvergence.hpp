// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>

#include <deque>

namespace sgpp {
namespace datadriven {

/**
 * A monitor to decide if a learning algorithm has converged. The
 * convergence criterion is based on the comparison of error
 * measurements throughout the training process.
 */
class RefinementMonitorConvergence : public RefinementMonitor {
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
  RefinementMonitorConvergence(double pDeclineThreshold, size_t pBufferSize,
                     size_t pMinRefInterval);
  /**
   * Destructor.
   */
  ~RefinementMonitorConvergence() override;

  /**
   * Stores the current error values in the buffer. If the buffer
   * has reached the maximum size, the oldest values are removed.
   *
   * @param numberInstances the number of instances that were used for the training step
   * @param currentValidError The current validation error
   * @param currentTrainError The current training error
   */
  void pushToBuffer(size_t numberInstances, double currentValidError,
      double currentTrainError) override;

  /**
   * Checks if the model needs to be refined. The convergence based monitor will at most
   * trigger one refinement at once.
   *
   * @return the number of refinements that are triggered by the monitor
   */
  size_t refinementsNecessary() override;

  /**
   * Examines the convergence criterion with the current
   * error observations.
   *
   * @return True if converged, false otherwise
   */

 private:
  bool checkConvergence();

  // counts how many iterations are yet to be performed until
  // next refinement can be triggered (only required if minRefInterval > 0 is
  // chosen)
  size_t nextRefCnt;
  size_t minRefInterval;
  // stores the latest validation error observations
  std::deque<double> validErrorDeclineBuffer;
  // stores the latest training error observations
  std::deque<double> trainErrorDeclineBuffer;
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

}  // namespace datadriven
}  // namespace sgpp

