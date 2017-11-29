// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>

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
 */
FloatArrayVector BSplineScalarProductEvaluator::get1DMixedIntegral(std::vector<double>& points,
                                                                   size_t index_i) {
  FloatArrayVector sums;
  for (size_t index_j = 0; index_j < points.size(); index_j++) {
    // performing Gauss-Legendre integration with twice as many points as for the simple integrals
    size_t numGaussPoints = 2 * ((degree + 1) / 2 + numAdditionalPoints);
    sgpp::base::DataVector roots;
    sgpp::base::DataVector quadratureweights;
    auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

    double sum = 0.0;
    std::vector<double> xi;
    createNakKnots(xValues, degree, xi);
    double productValue = 0.0;

    // constant function for single point, Lagrange polynomials for up to 9 points, B-splines
    // otherwise
    if (xValues.size() == 1) {
      numGaussPoints = 1;
      quadRule.getLevelPointsAndWeightsNormalized(
          std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
      sum = 1.0 * this->weight_function(roots[0]) * quadratureweights[0];
    } else if (xValues.size() < 9) {
      numGaussPoints = 2 * xValues.size();
      quadRule.getLevelPointsAndWeightsNormalized(
          std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x = roots[i];
        productValue =
            LagrangePolynomial(x, xValues, index_j) * LagrangePolynomial(x, xValues, index_i);
        double integrand = productValue * this->weight_function(x);
        sum += integrand * quadratureweights[i];
      }
    } else {
      quadRule.getLevelPointsAndWeightsNormalized(
          std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
      size_t first_segment_j = std::max(degree, index_j);
      size_t last_segment_j = std::min(xi.size() - degree - 1, index_j + degree + 1);
      size_t first_segment_i = std::max(degree, index_i);
      size_t last_segment_i = std::min(xi.size() - degree - 1, index_i + degree + 1);
      size_t first_segment = std::min(first_segment_i, first_segment_j);
      size_t last_segment = std::max(last_segment_i, last_segment_j);
      for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
        double a = std::max(0.0, xi[segmentIndex]);
        double b = std::min(1.0, xi[segmentIndex + 1]);
        double width = b - a;

        for (size_t i = 0; i < roots.getSize(); ++i) {
          double x = a + width * roots[i];
          productValue =
              nonUniformBSpline(x, degree, index_j, xi) * nonUniformBSpline(x, degree, index_i, xi);
          double integrand = productValue * this->weight_function(x);
          // multiply weights by length_old_interval / length_new_interval
          sum += integrand * quadratureweights[i] * width;
        }
      }
    }
    sums[index_j] = sum;
  }
  return sums;
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
void BSplineScalarProductEvaluator::calculate1DBSplineScalarProducts(
    std::vector<double>& points, std::vector<FloatArrayVector>& basisValues) {
  for (size_t index_i = 0; index_i < points.size(); ++index_i) {
    basisValues.push_back(get1DMixedIntegral(points, index_i));
  }
}

BSplineScalarProductEvaluator::~BSplineScalarProductEvaluator() {}

bool BSplineScalarProductEvaluator::needsOrderedPoints() { return true; }

bool BSplineScalarProductEvaluator::needsParameter() { return false; }

void BSplineScalarProductEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  basisValues.clear();
  calculate1DBSplineScalarProducts(xValues, basisValues);

  // is this ever used?
  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the basis values with the weight function
    for (size_t i = 0; i < basisValues.size(); ++i) {
      for (size_t j = 0; j < basisValues[i].size(); ++j) {
        // basisValues[i].scalarMult(weight_function(xValues[i]));
        sum += basisValues[i][j].getValue();
      }
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < basisValues.size(); ++i) {
      basisValues[i].scalarMult(sumInv);
    }
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector> >
BSplineScalarProductEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector> >(
      new BSplineScalarProductEvaluator(*this));
}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(3) {}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(degree) {}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints),
      degree(degree) {}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(
    BSplineScalarProductEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints),
      degree(other.degree) {}

void BSplineScalarProductEvaluator::setParameter(const FloatArrayVector& param) { return; }

void BSplineScalarProductEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
