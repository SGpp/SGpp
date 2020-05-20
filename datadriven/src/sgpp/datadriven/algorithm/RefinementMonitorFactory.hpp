// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Factory to create refinement monitors
 */
class RefinementMonitorFactory {
 public:
  RefinementMonitorFactory() = default;

  /**
   * Creates a refinement monitor
   * @param adaptivityConfig configuration for the adaptivity of the sparse grid
   * @return a new refinement monitor instance
   */
  RefinementMonitor* createRefinementMonitor(
      const sgpp::base::AdaptivityConfiguration& adaptivityConfig) const;
};

}  // namespace datadriven
}  // namespace sgpp
