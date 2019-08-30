// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/grid/distribution/ChebyshevDistribution.hpp>

#include <cmath>
#include <stdexcept>

namespace sgpp {
namespace combigrid {

ChebyshevDistribution::~ChebyshevDistribution() {}

double ChebyshevDistribution::compute(size_t numPoints, size_t j) {
  if (j >= numPoints) {
    throw std::logic_error("ChebyshevDistribution::compute: j >= numPoints");
  }

  return 0.5 * (1 - cos((static_cast<double>(2*j + 1) / static_cast<double>(2*numPoints)) * M_PI));
}

} /* namespace combigrid */
} /* namespace sgpp*/
