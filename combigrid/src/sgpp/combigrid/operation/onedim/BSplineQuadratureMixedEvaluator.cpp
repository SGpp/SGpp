// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureMixedEvaluator.hpp>
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
 * [ ToDo (rehmemk) but setting the segments depending on b_i AND b_j should speed up this
 * calculation!]
 * ToDo (rehmemk) simply added a loop over index_i to calculate FloatArrayVectors instead of
 * doubles. This vector-wise calculation can be optimized
 */
FloatArrayVector BSplineQuadratureMixedEvaluator::get1DMixedIntegral(std::vector<double>& points,
                                                                     size_t index_j) {
  FloatArrayVector sums;
  for (size_t index_i = 0; index_i < points.size(); index_i++) {
    // performing Gauss-Legendre integration with twice as many points as for the simple integrals
    size_t numGaussPoints = 2 * ((degree + 1) / 2 + numAdditionalPoints);
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
    sums[index_i] = sum;
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
void BSplineQuadratureMixedEvaluator::calculate1DMixedBSplineIntegrals(
    std::vector<double>& points, std::vector<FloatArrayVector>& basisValues) {
  size_t tempindex = 0;
  for (size_t index_j = 0; index_j < points.size(); ++index_j) {
    basisValues.push_back(get1DMixedIntegral(points, index_j));
    tempindex++;
  }
}

BSplineQuadratureMixedEvaluator::~BSplineQuadratureMixedEvaluator() {}

bool BSplineQuadratureMixedEvaluator::needsOrderedPoints() { return true; }

bool BSplineQuadratureMixedEvaluator::needsParameter() { return false; }

void BSplineQuadratureMixedEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  basisValues.clear();
  calculate1DMixedBSplineIntegrals(xValues, basisValues);

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
BSplineQuadratureMixedEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector> >(
      new BSplineQuadratureMixedEvaluator(*this));
}

BSplineQuadratureMixedEvaluator::BSplineQuadratureMixedEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(3) {}

BSplineQuadratureMixedEvaluator::BSplineQuadratureMixedEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      degree(degree) {}

BSplineQuadratureMixedEvaluator::BSplineQuadratureMixedEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints),
      degree(degree) {}

BSplineQuadratureMixedEvaluator::BSplineQuadratureMixedEvaluator(
    BSplineQuadratureMixedEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints),
      degree(other.degree) {}

void BSplineQuadratureMixedEvaluator::setParameter(const FloatArrayVector& param) { return; }

void BSplineQuadratureMixedEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
