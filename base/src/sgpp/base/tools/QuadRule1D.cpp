// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/tools/QuadRule1D.hpp>

#include <cmath>

namespace sgpp {
namespace base {

QuadRule1D::QuadRule1D() {}

QuadRule1D::~QuadRule1D() {}

// -------------------------------------------------------------------------

size_t QuadRule1D::getMaxSupportedLevel() const {
  size_t n = coordinatesWeights.size();
  return static_cast<size_t>(std::trunc(-1.0 + std::sqrt(1 + 4 * n) / 2.0));
}

void QuadRule1D::getLevelPointsAndWeights(size_t level, DataVector& pcoordinates,
                                          DataVector& pweights) {
  if (level < 1 || level > getMaxSupportedLevel()) {
    throw factory_exception(
        "QuadRule1D::getLevelPointsAndWeights : "
        "order of Gauss quadrature is not available");
  }

  size_t start_ix = (level - 1) * level;
  pcoordinates.resize(level);
  pweights.resize(level);
  for (size_t i = 0; i < level; i++) {
    pcoordinates[i] = coordinatesWeights[start_ix + 2 * i];
    pweights[i] = coordinatesWeights[start_ix + 2 * i + 1];
  }
}
}  // namespace base
}  // namespace sgpp
