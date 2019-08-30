/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * RefinementMonitorPeriodic.hpp
 *
 *  Created on: Jun 30, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>

#include <deque>

namespace sgpp {
namespace datadriven {

/**
 * A monitor that decides whether refinements should be performed using a simple periodic approach:
 * After at least a certain amount of instances has arrived, the monitor will allow a new refinement
 */
class RefinementMonitorPeriodic : public RefinementMonitor {
 public:
  /**
   * Constructor for the periodic refinement monitor
   *
   * @param period the number of instances that is needed to trigger a new refinement. If set to
   * zero, this monitor will always trigger exactly one refinement
   */
  explicit RefinementMonitorPeriodic(size_t period);
  /**
   * Destructor.
   */
  ~RefinementMonitorPeriodic() override;

  /**
   * Pushes a new iteration to the monitors state.
   *
   * @param numberInstances the number of instances that were used for the training step
   * @param currentValidError The current validation error
   * @param currentTrainError The current training error
   */
  void pushToBuffer(size_t numberInstances, double currentValidError,
      double currentTrainError) override;

  /**
   * Checks if the model needs to be refined. If multiples of the period data instances have arrived
   * this monitor might trigger more than one refinement
   *
   * @return the number of refinements that are triggered by the monitor
   */
  size_t refinementsNecessary() override;


 private:
  // Stores the number of instances required to trigger a new refinement
  size_t period;
  // Stores the total amount of instances that are currently being tracked
  size_t instancesTracked;
};

}  // namespace datadriven
}  // namespace sgpp

