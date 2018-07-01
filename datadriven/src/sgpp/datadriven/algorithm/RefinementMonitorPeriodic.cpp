/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * RefinementMonitorPeriodic.cpp
 *
 *  Created on: Jun 30, 2018
 *      Author: dominik
 */


#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>

namespace sgpp {
namespace datadriven {

RefinementMonitorPeriodic::RefinementMonitorPeriodic(size_t period)
    : period(0), instancesTracked(0) {}

RefinementMonitorPeriodic::~RefinementMonitorPeriodic() {}

void RefinementMonitorPeriodic::pushToBuffer(size_t numberInstances,
                                      double currentValidError,
                                      double currentTrainError) {
  (void)currentValidError;
  (void)currentTrainError;
  instancesTracked += numberInstances;
}

size_t RefinementMonitorPeriodic::refinementsNecessary() {
  // If there is no period given, this monitor will always allow refinements
  if (period == 0) return 1;

  size_t numRefinements = instancesTracked / period;
  instancesTracked -= numRefinements * period;
  return numRefinements;
}

}  // namespace datadriven
}  // namespace sgpp




