// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/utils/CombigridBSplineBasis.hpp>

#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator calculates the integrals int b_i(x) dx for B splines b_i. This is
 * done via quadrature based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class BSplineQuadratureEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  std::vector<double> xValues;
  std::vector<FloatScalarVector> basisValues;
  std::vector<double> basisCoefficients;
  sgpp::combigrid::SingleFunction weight_function;
  size_t numAdditionalPoints;
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t degree;
  double a;
  double b;

  /**
   * Calculates the integral of the B spline corresponding to the point with the given index
   * @param points grid points of the one dimensional grid
   * @param index index of the B-Spline whose integral will be calculated
   */
  double get1DIntegral(std::vector<double> &points, size_t index);

  /**
   * This function calculates the integral of the B spline basis given by their points
   * @param points The vector with the points
   * @param basisValues The integrals will be added to the back of this vector in the order of the
   * @param incrementQuadraturePoints increment for numAdditionalPoints in the iterative
   * quadrature
   * routine for custom weight function
   * @param tol tolerance for the iterative quadrature routine for custom weight function
   * points in the vector with the points. It is recommended to clear the basisValues vector before
   * calling this function to ensure that the basisValues are at the same position as their points
   */
  void calculate1DBSplineIntegrals(std::vector<double> &points,
                                   std::vector<FloatScalarVector> &basisValues,
                                   size_t incrementQuadraturePoints = 1, double tol = 1e-14);

 public:
  BSplineQuadratureEvaluator();
  explicit BSplineQuadratureEvaluator(size_t degree);

  /**
   * @param degree degree of the B splien basis
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
  BSplineQuadratureEvaluator(size_t degree, sgpp::combigrid::SingleFunction weight_function,
                             size_t numAdditionalPoints, bool normalizeWeights = true);

  BSplineQuadratureEvaluator(size_t degree, sgpp::combigrid::SingleFunction weight_function,
                             size_t numAdditionalPoints, double a, double b,
                             bool normalizeWeights = true);

  BSplineQuadratureEvaluator(BSplineQuadratureEvaluator const &other);
  ~BSplineQuadratureEvaluator() override;

  std::vector<FloatScalarVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;

  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;
  void setWeightFunction(sgpp::combigrid::SingleFunction weight_function) override {
    this->weight_function = weight_function;
    this->isCustomWeightFunction = true;
  }
  void setBounds(double a, double b) override {
    this->a = a;
    this->b = b;
  }

  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
};

} /* namespace combigrid */
} /* namespace sgpp*/
