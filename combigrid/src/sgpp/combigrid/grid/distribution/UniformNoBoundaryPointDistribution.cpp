// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "UniformNoBoundaryPointDistribution.hpp"

#include <stdexcept>

namespace sgpp {
namespace combigrid {

UniformNoBoundaryPointDistribution::~UniformNoBoundaryPointDistribution() {}

double UniformNoBoundaryPointDistribution::compute(size_t numPoints, size_t j) {
  if (j >= numPoints) {
    throw std::logic_error("UniformNoBoundaryPointDistribution::compute: j >= numPoints");
  }

  return static_cast<double>(j + 1) / static_cast<double>(numPoints + 1);
}

} /* namespace combigrid */
} /* namespace sgpp*/
