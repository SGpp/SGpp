// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
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

/**
 * Calculates the weight for the specific point
 */
double PolynomialQuadratureEvaluator::getWeight(std::vector<double>& points, size_t point) {
  LagrangePolynom p;
  p.points = points;
  p.point = point;
  size_t numGaussPoints = (p.degree() + 2) / 2 + numAdditionalPoints;

  return GaussLegendreQuadrature(numGaussPoints).evaluate([&p, this](double x) {
    return p.evaluate(x) * this->weight_function(x);
  });
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
void PolynomialQuadratureEvaluator::calculateWeights(std::vector<double>& points,
                                                     std::vector<FloatScalarVector>& weights) {
  // calc weight for each point
  for (size_t i = 0; i < points.size(); ++i) {
    weights.push_back(FloatScalarVector(getWeight(points, i)));
  }
}

PolynomialQuadratureEvaluator::~PolynomialQuadratureEvaluator() {}

bool PolynomialQuadratureEvaluator::needsOrderedPoints() { return false; }

bool PolynomialQuadratureEvaluator::needsParameter() { return false; }

void PolynomialQuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  weights.clear();
  calculateWeights(xValues, weights);

  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the weights with the weight function
    for (size_t i = 0; i < weights.size(); ++i) {
      // weights[i].scalarMult(weight_function(xValues[i]));
      sum += weights[i].getValue();
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(sumInv);
    }
  }
  /*else {
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i].scalarMult(weight_function(xValues[i]));
    }
  }*/
}

double PolynomialQuadratureEvaluator::getAbsoluteWeightSum() const {
  double abssum = 0.0;
  // double sum = 0.0;

  for (size_t i = 0; i < weights.size(); ++i) {
    // sum += weights[i].value();
    abssum += fabs(weights[i].value());
  }

  // std::cout << "sum: " << sum << "\n";

  return abssum;
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
PolynomialQuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new PolynomialQuadratureEvaluator(*this));
}

PolynomialQuadratureEvaluator::PolynomialQuadratureEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0) {}

PolynomialQuadratureEvaluator::PolynomialQuadratureEvaluator(
    sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints) {}

PolynomialQuadratureEvaluator::PolynomialQuadratureEvaluator(
    PolynomialQuadratureEvaluator const& other)
    : xValues(other.xValues),
      weights(other.weights),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints) {}

void PolynomialQuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

void PolynomialQuadratureEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

// CombiEvaluatorTypes PolynomialQuadratureEvaluator::getType() {
//  return CombiEvaluatorTypes::Scalar_PolynomialQuadrature;
//}

} /* namespace combigrid */
} /* namespace sgpp*/
