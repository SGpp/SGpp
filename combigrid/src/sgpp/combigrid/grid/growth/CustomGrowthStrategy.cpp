// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/growth/CustomGrowthStrategy.hpp>

namespace sgpp {
namespace combigrid {

CustomGrowthStrategy::CustomGrowthStrategy(std::function<size_t(size_t)> func) : func(func) {}

CustomGrowthStrategy::~CustomGrowthStrategy() {}

size_t CustomGrowthStrategy::numPoints(size_t level) { return func(level); }

} /* namespace combigrid */
} /* namespace sgpp*/
