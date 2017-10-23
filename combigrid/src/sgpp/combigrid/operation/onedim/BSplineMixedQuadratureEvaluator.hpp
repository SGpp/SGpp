// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSplineMIXEDQUADRATUREEVALUATOR_HPP_
#define BSplineMIXEDQUADRATUREEVALUATOR_HPP_

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator does quadrature based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class BSplineMixedQuadratureEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  std::vector<double> xValues;
  std::vector<FloatScalarVector> integrals;
  std::vector<double> basisCoefficients;
  sgpp::combigrid::SingleFunction weight_function;
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t numAdditionalPoints;  // additional gauss points used for a custom weight function
  size_t degree;

  double get1DMixedIntegral(std::vector<double> &points, size_t index_i, size_t index_j);

  void calculate1DMixedBSplineIntegrals(std::vector<double> &points,
                                        std::vector<FloatScalarVector> &integrals);

 public:
  BSplineMixedQuadratureEvaluator();
  BSplineMixedQuadratureEvaluator(size_t degree);

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
  BSplineMixedQuadratureEvaluator(size_t degree, sgpp::combigrid::SingleFunction weight_function,
                                  bool normalizeWeights = true, size_t numAdditionalPoints = 10);
  BSplineMixedQuadratureEvaluator(BSplineMixedQuadratureEvaluator const &other);
  virtual ~BSplineMixedQuadratureEvaluator();

  std::vector<FloatScalarVector> getBasisValues() override { return integrals; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setFunctionValuesAtGridPoints(std::vector<double> &functionValues) override;

  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;

  // The following is simply copied from QuadratureEvaluator. Applicable here?
  // can be used as a measure of stability of the quadrature algorithm. Minimum (and optimum) in
  // case of normalized weights is 1.0, i.e. all weights are non-negative.
  //  double getAbsoluteWeightSum() const;

  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
};

#endif /* BSplineMIXEDQUADRATUREEVALUATOR_HPP_ */

} /* namespace combigrid */
} /* namespace sgpp*/
