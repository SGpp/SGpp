// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>

#include <iomanip>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Struct for a LagrangePolynom, used to eval it
 */
struct LagrangePolynom {
  std::vector<double> points;
  size_t point;

  double evaluate(double x) {
    double result = 1.0;
    for (size_t i = 0; i < points.size(); ++i) {
      if (i != point) {
        // TODO(holzmudd): precalculate denominator?
        result *= (x - points[i]) / (points[point] - points[i]);
      }
    }
    return result;
  }

  size_t degree() { return points.size(); }
};

double doGaussLegendreQuadrature(double a, double b, LagrangePolynom& f) {
  // load gauss weights and roots for Gauss Legendre quadrature
  size_t level = f.degree() + 1;
  base::DataVector roots;
  base::DataVector weights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  quadRule.getLevelPointsAndWeightsNormalized(level, roots, weights);

  // evaluate the Lagrange polynomial at the roots of the
  // Legendre polynomials
  double width = b - a;
  double sum = 0.0;
  for (size_t i = 0; i < roots.getSize(); ++i) {
    sum += weights[i] * f.evaluate(width * roots[i] + a);
  }
  return width * sum;
}

/**
 * Integrates the polynom
 */
double integrate(LagrangePolynom polynom) {
  // return romberg2(polynom);
  return doGaussLegendreQuadrature(0.0, 1.0, polynom);
}

/**
 * Calculates the weight for the specific point
 */
double getWeight(std::vector<double>& points, size_t point) {
  LagrangePolynom p;
  p.points = points;
  p.point = point;

  return integrate(p);
}

/**
 * This Function calculates the weights of the given points, each weight is calculated individually
 * @param points The vector with the points, they dont need to have a specific order
 * @param weights The weights will be added to the back of this vector in the order of the points in
 * the vector with the points,
 * it is recommended to clear the weight vector before calling this function to ensure that the
 * weights are at the same position
 * as their points
 */
void calculateWeights(std::vector<double>& points, std::vector<FloatScalarVector>& weights) {
  // calc weight for each point
  for (size_t i = 0; i < points.size(); ++i) {
    weights.push_back(FloatScalarVector(getWeight(points, i)));
  }
}

QuadratureEvaluator::~QuadratureEvaluator() {}

bool QuadratureEvaluator::needsOrderedPoints() { return false; }

bool QuadratureEvaluator::needsParameter() { return false; }

void QuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  weights.clear();
  calculateWeights(xValues, weights);

  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the weights with the weight function
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(weight_function(xValues[i]));
      sum += weights[i].getValue();
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(sumInv);
    }
  } else {
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(weight_function(xValues[i]));
    }
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> > QuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new QuadratureEvaluator(*this));
}

QuadratureEvaluator::QuadratureEvaluator(sgpp::combigrid::SingleFunction weight_function,
                                         bool normalizeWeights)
    : weight_function(weight_function), normalizeWeights(normalizeWeights) {}

QuadratureEvaluator::QuadratureEvaluator(QuadratureEvaluator const& other)
    : xValues(other.xValues),
      weights(other.weights),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights) {}

void QuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

} /* namespace combigrid */
} /* namespace sgpp*/
