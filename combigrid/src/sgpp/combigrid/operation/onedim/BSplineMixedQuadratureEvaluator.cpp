// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineMixedQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Calculates the weight for the specific point
 * @param points grid points of the one dimensional grid the interpolation will be performed on
 * @param index_i index of B-spline b_i
 * @param index_j index of B-spline b_j
 * @return integral of b_i*b_j
 *
 * w.l.o.g. integrate over support of b_i because b_i*b_j is zero outside supp(b_i)
 * [ but setting the segments depending on b_i AND b_j should speed up this calculation!]
 */
double BSplineMixedQuadratureEvaluator::get1DMixedIntegral(std::vector<double>& points,
                                                           size_t index_i, size_t index_j) {
  // performing Gauss-Legendre integration
  size_t numGaussPoints = (degree + 1) / 2 + numAdditionalPoints;
  sgpp::base::DataVector roots;
  sgpp::base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  double sum = 0.0;
  std::vector<double> xi;
  createNakKnots(xValues, degree, xi);
  double bsplinevalue = 0.0;

  if (xValues.size() == 1) {
    numGaussPoints = 1;
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    sum = 1.0 * this->weight_function(roots[0]) * quadratureweights[0];
  } else if (xValues.size() < 9) {
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x = roots[i];
      bsplinevalue =
          LagrangePolynomial(x, xValues, index_i) * LagrangePolynomial(x, xValues, index_j);
      double integrand = bsplinevalue * this->weight_function(x);
      sum += integrand * quadratureweights[i];
    }
  } else {
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    // use segments of b_i. Ass bi*bj =0 outside supp(b_i) this works fine.
    // ToDo (rehmemk) customizing first-/last_segment to bi AND bj would lead to a speed up!
    size_t first_segment = std::max(degree, index_i);
    size_t last_segment = std::min(xi.size() - degree - 1, index_i + degree + 1);
    for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
      double a = std::max(0.0, xi[segmentIndex]);
      double b = std::min(1.0, xi[segmentIndex + 1]);
      double width = b - a;

      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x = a + width * roots[i];
        bsplinevalue =
            nonUniformBSpline(x, degree, index_i, xi) * nonUniformBSpline(x, degree, index_j, xi);
        double integrand = bsplinevalue * this->weight_function(x);
        // multiply weights by length_old_interval / length_new_interval
        sum += integrand * quadratureweights[i] * width;
      }
    }
  }
  return sum;
}

/**
 * This Function calculates the weights of the given points, each weight is calculated
 * individually
 * @param points The vector with the points, they dont need to have a specific order
 * @param integrals The integrals will be added to the back of this vector in the order of the
 * points in the vector with the points,
 * it is recommended to clear the weight vector before calling this function to ensure that the
 * weights are at the same position
 * as their points
 */
void BSplineMixedQuadratureEvaluator::calculate1DMixedBSplineIntegrals(
    std::vector<double>& points, std::vector<FloatScalarVector>& integrals) {
  // "weights" here are the integrals!
  std::cout << "integrals:\n";
  size_t tempindex = 0;
  for (size_t index_i = 0; index_i < points.size(); ++index_i) {
    for (size_t index_j = 0; index_j < points.size(); ++index_j) {
      integrals.push_back(FloatScalarVector(get1DMixedIntegral(points, index_i, index_j)));
      std::cout << integrals[tempindex] << " ";
      tempindex++;
    }
  }
  std::cout << "\n";
}

BSplineMixedQuadratureEvaluator::~BSplineMixedQuadratureEvaluator() {}

bool BSplineMixedQuadratureEvaluator::needsOrderedPoints() { return true; }

bool BSplineMixedQuadratureEvaluator::needsParameter() { return false; }

void BSplineMixedQuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  integrals.clear();
  calculate1DMixedBSplineIntegrals(xValues, integrals);

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

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
BSplineMixedQuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new BSplineMixedQuadratureEvaluator(*this));
}

BSplineMixedQuadratureEvaluator::BSplineMixedQuadratureEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(3) {}

BSplineMixedQuadratureEvaluator::BSplineMixedQuadratureEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(degree) {}

BSplineMixedQuadratureEvaluator::BSplineMixedQuadratureEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints),
      degree(degree) {}

BSplineMixedQuadratureEvaluator::BSplineMixedQuadratureEvaluator(
    BSplineMixedQuadratureEvaluator const& other)
    : xValues(other.xValues),
      integrals(other.integrals),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints),
      degree(other.degree) {}

void BSplineMixedQuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

void BSplineMixedQuadratureEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
