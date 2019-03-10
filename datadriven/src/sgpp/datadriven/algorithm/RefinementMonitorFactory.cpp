/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RefinementMonitorFactory.cpp
 *
 *  Created on: Jul 1, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>

namespace sgpp {
namespace datadriven {

RefinementMonitor* RefinementMonitorFactory::createRefinementMonitor(
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig) const {
  if (adaptivityConfig.errorBasedRefinement) {
    return new RefinementMonitorConvergence(adaptivityConfig.errorConvergenceThreshold,
        adaptivityConfig.errorBufferSize, adaptivityConfig.errorMinInterval);
  } else {
    return new RefinementMonitorPeriodic(adaptivityConfig.refinementPeriod);
  }
}
//RefinementMonitor* RefinementMonitorFactory::createCoarseningMonitor(
//            const sgpp::base::CoarseningConfiguration& coarseningConfig) const {
//  if (coarseningConfig.errorBasedCoarsening) {
//    return new RefinementMonitorConvergence(coarseningConfig.errorConvergenceThreshold,
//                                            coarseningConfig.errorBufferSize, coarseningConfig.errorMinInterval);
//  } else {
//    return new RefinementMonitorPeriodic(coarseningConfig.CoarseningPeriod);
//  }
//}
}  // namespace datadriven
}  // namespace sgpp

