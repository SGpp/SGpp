// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "UniformPointDistribution.hpp"

#include <stdexcept>

namespace sgpp {
namespace combigrid {

UniformPointDistribution::~UniformPointDistribution() {}

double UniformPointDistribution::compute(size_t numPoints, size_t j) {
  if (j >= numPoints) {
    throw std::logic_error("UniformPointDistribution::compute: j >= numPoints");
  }

  if (numPoints == 1) {
    return 0.5;
  }

  return static_cast<double>(j) / static_cast<double>(numPoints - 1);
}

} /* namespace combigrid */
} /* namespace sgpp*/
