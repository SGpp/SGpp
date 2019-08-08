// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp>

#include <cmath>
#include <stdexcept>

namespace sgpp {
namespace combigrid {

ClenshawCurtisDistribution::~ClenshawCurtisDistribution() {}

double ClenshawCurtisDistribution::compute(size_t numPoints, size_t j) {
  if (j >= numPoints) {
    throw std::logic_error("ClenshawCurtisDistribution::compute: j >= numPoints");
  }

  if (numPoints == 1) {
    return 0.5;
  }

  return 0.5 * (1 - cos((static_cast<double>(j) / static_cast<double>(numPoints - 1)) * M_PI));
}

} /* namespace combigrid */
} /* namespace sgpp*/
