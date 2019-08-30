// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSHERMITEQUADRULE1D_HPP_
#define GAUSSHERMITEQUADRULE1D_HPP_

#include <sgpp/base/tools/QuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * load gauss quadrature points for standard normal weight function. The points
 * and the weights are generated with numpy.polynomial.hermite.hermgauss,
 * the coordinates are scaled by sqrt(2), the weights are normalized to 1.
 */
class GaussHermiteQuadRule1D : public QuadRule1D {
 public:
  GaussHermiteQuadRule1D();
  ~GaussHermiteQuadRule1D() override;

  // delete the copy constructor
  GaussHermiteQuadRule1D(const GaussHermiteQuadRule1D& that) = delete;

  /**
   * the coordinates are scaled by sqrt(2) and then normalized with respect
   * to a given mean and standard deviation. The weights are normalized
   * to 1.
   *
   * @param level level of quadrature, is equal to the number of quadrature points
   * @param coordinates returns the x-coordinates in [-infty, infty]
   * @param weights returns the corresponding weights (scaled by sqrt(2))
   * @param mean mean of the normal distribution the coordinates should be transformed to
   * @param stdd standard deviation of the normal distribution the coordinates should be transformed
   * to
   */
  void getLevelPointsAndWeightsNormalized(size_t level, base::DataVector& coordinates,
                                          base::DataVector& weights, double mean = 0.0,
                                          double stdd = 1.0);

  static GaussHermiteQuadRule1D& getInstance();
};

}  // namespace base
}  // namespace sgpp

#endif /* GAUSSHERMITEQUADRULE1D_HPP_ */
