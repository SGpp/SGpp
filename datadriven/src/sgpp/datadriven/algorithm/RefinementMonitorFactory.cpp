// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>

namespace sgpp {
namespace datadriven {

RefinementMonitor* RefinementMonitorFactory::createRefinementMonitor(
    const sgpp::base::AdaptivityConfiguration& adaptConfig) const {
  if (adaptConfig.errorBasedRefinement_) {
    return new RefinementMonitorConvergence(adaptConfig.errorConvergenceThreshold_,
        adaptConfig.errorBufferSize_, adaptConfig.errorMinInterval_);
  } else {
    return new RefinementMonitorPeriodic(adaptConfig.refinementPeriod_);
  }
}

}  // namespace datadriven
}  // namespace sgpp

