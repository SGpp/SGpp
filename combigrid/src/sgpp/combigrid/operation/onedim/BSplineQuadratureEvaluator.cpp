// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Calculates the weight for the specific point
 * @param points grid points of theone dimensional grid the interpolation will be performed on
 * @param index index of the B-Spline whose integral will be calculated
 */
double BSplineQuadratureEvaluator::getIntegral(std::vector<double>& points, size_t index) {
  // performing Gauss-Legendre integration
  size_t numGaussPoints = (degree + 1) / 2 + numAdditionalPoints;
  base::DataVector roots;
  base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  quadRule.getLevelPointsAndWeightsNormalized(
      std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);

  double sum = 0.0;
  std::vector<double> xi;
  createNakKnots(xValues, degree, xi);

  // multiply weights by length_old_interval / length_new_interval
  for (size_t i = 0; i < quadratureweights.size(); i++) {
    quadratureweights[i] *= 1.0 / (points.size() - 1.0);
  }

  size_t offset = (degree + 1) / 2;
  size_t first_segment = std::max(offset, index);
  size_t last_segment = std::min(xi.size() - offset - 1, index + degree + 1);
  for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
    double a = std::max(0.0, xi[segmentIndex]);
    double b = std::min(1.0, xi[segmentIndex + 1]);
    double width = b - a;

    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x = a + width * roots[i];
      double bsplinevalue = nonUniformBSpline(x, degree, index, xi);

      double integrand = bsplinevalue * this->weight_function(x);
      sum += integrand * quadratureweights[i];
    }
  }
  return sum;
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
void BSplineQuadratureEvaluator::calculateBSplineIntegrals(
    std::vector<double>& points, std::vector<FloatScalarVector>& integrals) {
  // "weights" here are the integrals!
  for (size_t index = 0; index < points.size(); ++index) {
    integrals.push_back(FloatScalarVector(getIntegral(points, index)));
  }
}

BSplineQuadratureEvaluator::~BSplineQuadratureEvaluator() {}

bool BSplineQuadratureEvaluator::needsOrderedPoints() { return true; }

bool BSplineQuadratureEvaluator::needsParameter() { return false; }

void BSplineQuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  integrals.clear();
  calculateBSplineIntegrals(xValues, integrals);

  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the weights with the weight function
    for (size_t i = 0; i < integrals.size(); ++i) {
      // weights[i].scalarMult(weight_function(xValues[i]));
      sum += integrals[i].getValue();
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < integrals.size(); ++i) {
      integrals[i].scalarMult(sumInv);
    }
  }
}
// The following is simply copied from QuadratureEvaluator. Applicable here?
// double BSplineQuadratureEvaluator::getAbsoluteWeightSum() const {
//  double abssum = 0.0;
//  // double sum = 0.0;
//
//  for (size_t i = 0; i < weights.size(); ++i) {
//    // sum += weights[i].value();
//    abssum += fabs(weights[i].value());
//  }
//
//  // std::cout << "sum: " << sum << "\n";
//
//  return abssum;
//}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
BSplineQuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new BSplineQuadratureEvaluator(*this));
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(3) {}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(degree) {}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints),
      degree(degree) {}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(BSplineQuadratureEvaluator const& other)
    : xValues(other.xValues),
      integrals(other.integrals),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints),
      degree(other.degree) {}

void BSplineQuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

void BSplineQuadratureEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
