// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <functional>
#include <vector>
#include <map>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator calculates the scalar products \f$\int b_i b_j\f$. This is done via quadrature
 * based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class PolynomialScalarProductEvaluator : public AbstractLinearEvaluator<FloatArrayVector> {
  std::vector<double> xValues;
  std::vector<FloatArrayVector> basisValues;
  std::vector<double> basisCoefficients;
  sgpp::combigrid::SingleFunction weight_function;
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t numAdditionalPoints;  // additional gauss points used for a custom weight function

  double xlower;
  double xupper;

  // lookup table for inner products
  std::map<size_t, double> scalarProductsMap;

  /**
   * Generate an unique key for two basis functions for reusing the scalar product
   *
   * @param i index of polynomial
   * @param j index of polynomial
   * @return unique key for hash map
   */
  size_t generateKey(size_t i, size_t j);

  /**
   * Performs Gauss-Legendre quadrature for the product of the given polynomials and the
   * weight function
   *
   * @param p_i Lagrange polynomial for point i
   * @param p_j Lagrange polynomial for point j
   * @return \f$\int p_i(x) p_j(x) f(x) dx\f$
   */
  double quad(LagrangePolynom &p_i, LagrangePolynom &p_j);

  /**
   * Calculates the weight for the specific point
   * @param points grid points of the one dimensional grid the interpolation will be performed on
   * @param index_i index of polynomial
   * @return integral of \f$b_i*b_j\f$
   */
  FloatArrayVector get1DMixedIntegral(std::vector<double> &points, size_t index_i);

  /**
   * This Function calculates the weights of the given points, each weight is calculated
   * individually
   * @param points The vector with the points, they don't need to have a specific order
   * @param integrals The integrals will be added to the back of this vector in the order of the
   * points in the vector with the points,
   * it is recommended to clear the weight vector before calling this function to ensure that the
   * weights are at the same position
   * as their points
   */
  void calculate1DPolynomialScalarProducts(std::vector<double> &points,
                                           std::vector<FloatArrayVector> &integrals);

 public:
  PolynomialScalarProductEvaluator();

  /**
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
  PolynomialScalarProductEvaluator(sgpp::combigrid::SingleFunction weight_function,
                                   bool normalizeWeights = true, size_t numAdditionalPoints = 10);
  PolynomialScalarProductEvaluator(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> orthogBasis);

  PolynomialScalarProductEvaluator(PolynomialScalarProductEvaluator const &other);
  ~PolynomialScalarProductEvaluator() override;

  std::vector<FloatArrayVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;

  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatArrayVector const &param) override;

  void initializeBounds(std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> orthogBasis);

  // The following is simply copied from QuadratureEvaluator. Applicable here?
  // can be used as a measure of stability of the quadrature algorithm. Minimum (and optimum) in
  // case of normalized weights is 1.0, i.e. all weights are non-negative.
  //  double getAbsoluteWeightSum() const;

  std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> cloneLinear() override;
};

} /* namespace combigrid */
} /* namespace sgpp*/
