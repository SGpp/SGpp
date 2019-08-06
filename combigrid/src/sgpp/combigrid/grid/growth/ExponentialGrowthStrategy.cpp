// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp>

namespace sgpp {
namespace combigrid {

ExponentialGrowthStrategy::~ExponentialGrowthStrategy() {}

size_t ExponentialGrowthStrategy::numPoints(size_t level) {
  return level == 0 ? 1 : (static_cast<size_t>(1) << level) + 1;
}

} /* namespace combigrid */
} /* namespace sgpp*/
