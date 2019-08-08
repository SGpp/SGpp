// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>

namespace sgpp {
namespace combigrid {

LinearGrowthStrategy::LinearGrowthStrategy(size_t factor) : factor(factor) {}

LinearGrowthStrategy::~LinearGrowthStrategy() {}

size_t LinearGrowthStrategy::numPoints(size_t level) { return 1 + factor * level; }

} /* namespace combigrid */
} /* namespace sgpp*/
