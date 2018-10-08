// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>

namespace sgpp {
namespace optimization {
double quad(std::function<double(double)> f, double a, double b, size_t quadOrder) {
  base::DataVector coordinates, weights;
  base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  // scale coordinates from [0,1] to [a,b]
  coordinates.mult(b - a);
  base::DataVector aVector(coordinates.getSize(), a);
  coordinates.add(aVector);

  double res = 0;
  for (size_t i = 0; i < coordinates.getSize(); i++) {
    res += weights[i] * f(coordinates[i]);
  }
  return res * (b - a);
}

}  // namespace optimization
}  // namespace sgpp
