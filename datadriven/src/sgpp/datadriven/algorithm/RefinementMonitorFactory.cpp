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
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig) const {
  if (adaptivityConfig.errorBasedRefinement_) {
    return new RefinementMonitorConvergence(adaptivityConfig.errorConvergenceThreshold_,
                                            adaptivityConfig.errorBufferSize_,
                                            adaptivityConfig.errorMinInterval_);
  } else {
    return new RefinementMonitorPeriodic(adaptivityConfig.refinementPeriod_);
  }
}

}  // namespace datadriven
}  // namespace sgpp
