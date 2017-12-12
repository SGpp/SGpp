// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator calculates the scalar products int b_i(x) b_j(x) dx for B splines b_i and b_j.
 * This is done via quadrature based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class BSplineScalarProductEvaluator : public AbstractLinearEvaluator<FloatArrayVector> {
  std::vector<double> xValues;
  std::vector<FloatArrayVector> basisValues;
  std::vector<double> basisCoefficients;
  sgpp::combigrid::SingleFunction weight_function;
  size_t numAdditionalPoints;  // additional gauss points used for a custom weight function
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t degree;

  /**
   * Calculates the one dimensional integrals int b_i(x) b_j(x) dx  for all j
   * @param points grid points of the one dimensional grid
   * @param index_i index of B-spline b_i
   * @return integral of b_i*b_j
   */
  FloatArrayVector get1DL2ScalarProduct(std::vector<double> const &points, size_t const &index_i);

  /**
   * This Function calculates all integrals int b_i(x) b_j(x) dx
   * @param points grid points of the one dimensional grid
   * @param integrals The integrals will be added to the back of this vector in the order of the
   * points in the vector with the points,
   * it is recommended to clear the weight vector before calling this function to ensure that the
   * weights are at the same position as their points
   */
  void calculate1DBSplineScalarProducts(std::vector<double> &points,
                                        std::vector<FloatArrayVector> &integrals);

 public:
  BSplineScalarProductEvaluator();
  explicit BSplineScalarProductEvaluator(size_t degree);

  /**
   * @param degree degree of the B spline basis
   * @param numAdditionalPoints Specifies how many Gauss-Legrendre points should be used in addition
   * to the default when integrating the Lagrange polynomials for computing the quadrature weights.
   * This number should be higher if the weight function is hard to integrate.
   * @param weight_function optional weight function w that can be included in the quadrature
   * weights. The Quadrature evaluator then approximates the integral of f*w instead of the integral
   * of f. This provides more precision by calling w more often than f, which might be useful if w
   * can be evaluated much faster than f.
   * @param normalizeWeights If this is set to true, the weights are normalized so that they sum up
   * to 1. This might be useful if the weight function is (or should be) a probability distribution
   * on the domain.
   */
  BSplineScalarProductEvaluator(size_t degree, sgpp::combigrid::SingleFunction weight_function,
                                size_t numAdditionalPoints, bool normalizeWeights = true);
  BSplineScalarProductEvaluator(BSplineScalarProductEvaluator const &other);
  virtual ~BSplineScalarProductEvaluator();

  std::vector<FloatArrayVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;

  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatArrayVector const &param) override;

  std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> cloneLinear() override;
};

} /* namespace combigrid */
} /* namespace sgpp*/
