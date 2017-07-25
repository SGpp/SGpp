// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRATUREEVALUATOR_HPP_
#define QUADRATUREEVALUATOR_HPP_

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
class QuadratureEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  std::vector<double> xValues;
  std::vector<FloatScalarVector> weights;
  sgpp::combigrid::SingleFunction weight_function;
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t numAdditionalPoints;  // additional gauss points used for a custom weight function

  double getWeight(std::vector<double> &points, size_t point);
  void calculateWeights(std::vector<double> &points, std::vector<FloatScalarVector> &weights);

 public:
  QuadratureEvaluator();
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
  QuadratureEvaluator(sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights = true,
                      size_t numAdditionalPoints = 10);
  QuadratureEvaluator(QuadratureEvaluator const &other);
  virtual ~QuadratureEvaluator();

  virtual std::vector<FloatScalarVector> getBasisCoefficients() { return weights; }

  virtual void setGridPoints(std::vector<double> const &newXValues);

  virtual bool needsOrderedPoints();
  virtual bool needsParameter();
  virtual void setParameter(FloatScalarVector const &param);

  // can be used as a measure of stability of the quadrature algorithm. Minimum (and optimum) in
  // case of normalized weights is 1.0, i.e. all weights are non-negative.
  double getAbsoluteWeightSum() const;

  virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
};

#endif /* QUADRATUREEVALUATOR_HPP_ */

} /* namespace combigrid */
} /* namespace sgpp*/
