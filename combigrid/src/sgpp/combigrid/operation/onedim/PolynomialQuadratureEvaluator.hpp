// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Class for a LagrangePolynom, used to eval it
 */
class LagrangePolynom {
 private:
  std::vector<double> points;
  size_t point;
  std::vector<double> denominator;

 public:
  void initialize(size_t newPoint, std::vector<double> &newPoints) {
    points = newPoints;
    point = newPoint;
    denominator.resize(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
      if (i != point) {
        denominator[i] = 1. / (points[point] - points[i]);
      } else {
        denominator[i] = 0.0;
      }
    }
  }

  double evaluate(double x) {
    double result = 1.0;
    for (size_t i = 0; i < points.size(); ++i) {
      if (i != point) {
        result *= (x - points[i]) * denominator[i];
      }
    }
    return result;
  }

  size_t degree() { return points.size(); }
};

/**
 * This evaluator does quadrature based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class PolynomialQuadratureEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  std::vector<double> xValues;
  std::vector<FloatScalarVector> weights;
  std::vector<double> basisCoefficients;
  sgpp::combigrid::SingleFunction weight_function;
  bool normalizeWeights;
  bool isCustomWeightFunction;
  size_t numAdditionalPoints;  // additional gauss points used for a custom weight function

  double getWeight(std::vector<double> &points, size_t point, GaussLegendreQuadrature &quadRule);
  void calculateWeights(std::vector<double> &points, std::vector<FloatScalarVector> &weights);

 public:
  PolynomialQuadratureEvaluator();
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
  PolynomialQuadratureEvaluator(sgpp::combigrid::SingleFunction weight_function,
                                bool normalizeWeights = true, size_t numAdditionalPoints = 10);
  PolynomialQuadratureEvaluator(PolynomialQuadratureEvaluator const &other);
  ~PolynomialQuadratureEvaluator() override;

  std::vector<FloatScalarVector> getBasisValues() override { return weights; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;

  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;

  // can be used as a measure of stability of the quadrature algorithm. Minimum (and optimum) in
  // case of normalized weights is 1.0, i.e. all weights are non-negative.
  double getAbsoluteWeightSum() const;
  double getRelativeConditionNumber() const;

  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
};

} /* namespace combigrid */
} /* namespace sgpp*/
