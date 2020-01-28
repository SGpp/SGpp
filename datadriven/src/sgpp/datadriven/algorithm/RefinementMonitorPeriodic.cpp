// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

RefinementMonitorPeriodic::RefinementMonitorPeriodic(size_t period)
    : period(period), instancesTracked(0) {}

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




