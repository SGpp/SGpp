// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Uses Gauss-Legendre quadrature to approximate the integral of a function on the domain [0, 1].
 */
class GaussLegendreQuadrature {
  base::DataVector roots;
  base::DataVector weights;

 public:
  /**
   * Constructor.
   * @param numPoints The number of grid points that should be used. This number should be > 0. If
   * not enough grid points are available, the maximum number of available grid points is taken.
   */
  explicit GaussLegendreQuadrature(size_t numPoints);

  void initialize(size_t numPoints);

  template <typename Func>
  double evaluate(Func const &func, double a = 0.0, double b = 1.0) {
    double width = b - a;
    double sum = 0.0;
    for (size_t i = 0; i < roots.getSize(); ++i) {
      sum += weights[i] * func(a + width * roots[i]);
    }
    return width * sum;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
