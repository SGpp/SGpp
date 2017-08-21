// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/tools/QuadRule1D.hpp>

namespace sgpp {
namespace base {

QuadRule1D::QuadRule1D(size_t maxSupportedLevel)
    : maxSupportedLevel(maxSupportedLevel),
      coordinates(maxSupportedLevel),
      weights(maxSupportedLevel) {}

QuadRule1D::~QuadRule1D() {
  for (size_t i = 0; i < coordinates.size(); i++) {
    delete coordinates[i];
  }

  for (size_t i = 0; i < weights.size(); i++) {
    delete weights[i];
  }
}

// -------------------------------------------------------------------------

void QuadRule1D::getLevelPointsAndWeights(size_t level, DataVector& pcoordinates,
                                          DataVector& pweights) {
  if (level < 1 || level > maxSupportedLevel) {
    throw factory_exception(
        "QuadRule1D::getLevelPointsAndWeights : "
        "order of gauss quadrature has to be within {1, ..., 50}");
  }

  pcoordinates = *coordinates[level - 1];
  pweights = *weights[level - 1];
}

size_t QuadRule1D::getMaxSupportedLevel() const { return maxSupportedLevel; }

}  // namespace base
}  // namespace sgpp
