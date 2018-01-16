// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

double BSplineQuadratureEvaluator::get1DIntegral(std::vector<double>& points, size_t index) {
  // performing Gauss-Legendre integration. Polynomials of degree 2*numGaussPoints-1 are integrated
  // exact
  size_t numGaussPoints = (degree + 1) / 2 + numAdditionalPoints;
  base::DataVector roots;
  base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();

  double sum = 0.0;
  std::vector<double> xi = createNakKnots(xValues, degree);
  double bsplinevalue = 0.0;

  // constant function for single point, Lagrange polynomials while not enough knots for not a
  // knot B-splines, nak B-splines otherwise
  if (xValues.size() == 1) {
    numGaussPoints = 1;
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    sum = 1.0 * this->weight_function(roots[0]) * quadratureweights[0];
  } else if ((degree == 3 && (xValues.size() < 5)) || ((degree == 5) && (xValues.size() < 9))) {
    numGaussPoints = xValues.size();
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x = roots[i];
      bsplinevalue = LagrangePolynomial(x, xValues, index);
      double integrand = bsplinevalue * this->weight_function(x);
      sum += integrand * quadratureweights[i];
    }
  } else {
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    size_t first_segment = std::max(degree, index);
    size_t last_segment = std::min(xValues.size(), index + degree + 1);
    for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
      double a = std::max(0.0, xi[segmentIndex]);
      double b = std::min(1.0, xi[segmentIndex + 1]);
      double width = b - a;

      // ToDo(rehmemk) Use GaussLegendreQuadrature.cpp's evaluate_iteratively if
      // weight function != 1
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x = a + width * roots[i];
        // ToDO(rehmemk) Rewrite this whole  routine , don't use createKnots and use the
        // Lagrange polynomials inside expUuniformNakBspline
        bsplinevalue = expUniformNakBspline(x, degree, index, xValues);
        double integrand = bsplinevalue * this->weight_function(x);
        // multiply weights by length_old_interval / length_new_interval
        sum += integrand * quadratureweights[i] * width;
      }
    }
  }
  return sum;
}

void BSplineQuadratureEvaluator::calculate1DBSplineIntegrals(
    std::vector<double>& points, std::vector<FloatScalarVector>& basisValues) {
  basisValues.resize(points.size());
  // #pragma omp parallel for schedule(static)
  for (size_t index = 0; index < points.size(); ++index) {
    basisValues[index] = FloatScalarVector(get1DIntegral(points, index));
  }
}

BSplineQuadratureEvaluator::~BSplineQuadratureEvaluator() {}

bool BSplineQuadratureEvaluator::needsOrderedPoints() { return true; }

bool BSplineQuadratureEvaluator::needsParameter() { return false; }

void BSplineQuadratureEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  basisValues.clear();
  calculate1DBSplineIntegrals(xValues, basisValues);

  if (normalizeWeights) {
    double sum = 0.0;

    // multiply the weights with the weight function
    for (size_t i = 0; i < basisValues.size(); ++i) {
      // weights[i].scalarMult(weight_function(xValues[i]));
      sum += basisValues[i].getValue();
    }

    double sumInv = 1.0 / sum;

    for (size_t i = 0; i < basisValues.size(); ++i) {
      basisValues[i].scalarMult(sumInv);
    }
  }
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
BSplineQuadratureEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new BSplineQuadratureEvaluator(*this));
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      numAdditionalPoints(0),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      degree(3) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = 3;
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      numAdditionalPoints(0),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      degree(degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = degree;
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
    bool normalizeWeights)
    : weight_function(weight_function),
      numAdditionalPoints(numAdditionalPoints),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      degree(degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = degree;
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(BSplineQuadratureEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      numAdditionalPoints(other.numAdditionalPoints),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      degree(other.degree) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = other.degree;
}

void BSplineQuadratureEvaluator::setParameter(const FloatScalarVector& param) { return; }

void BSplineQuadratureEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
