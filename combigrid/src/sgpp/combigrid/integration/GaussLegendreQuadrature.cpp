// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

#include <algorithm>

namespace sgpp {
namespace combigrid {

GaussLegendreQuadrature::GaussLegendreQuadrature(size_t numPoints) { initialize(numPoints); }

void GaussLegendreQuadrature::initialize(size_t numPoints) {
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  if (numPoints > 0) {
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numPoints, quadRule.getMaxSupportedLevel()), roots, weights);
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
