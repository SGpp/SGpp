// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

namespace sgpp {
namespace datadriven {

/**
 * Superclass for refinement monitors. They track whether a refinement should happen or not.
 */
class RefinementMonitor {
 public:
  /**
   * Destructor
   */
  virtual ~RefinementMonitor() = default;

  /**
   * Stores the current error values in the buffer. If the buffer
   * has reached the maximum size, the oldest values are removed.
   *
   * @param numberInstances the number of instances added
   * @param currentValidError The current validation error
   * @param currentTrainError The current training error
   */
  virtual void pushToBuffer(size_t numberInstances, double currentValidError,
      double currentTrainError) = 0;
  /**
   * Checks if the model needs to be refined.
   *
   * @return the number of refinements that are triggered by the monitor
   */
  virtual size_t refinementsNecessary() = 0;
};

}  // namespace datadriven
}  // namespace sgpp

