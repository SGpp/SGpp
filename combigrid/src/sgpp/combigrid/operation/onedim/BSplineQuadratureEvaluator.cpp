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
  //  std::cout << "BsplineQuadratureEvaluator: " << std::endl;
  //  std::cout << "index " << index << std::endl;
  // performing Gauss-Legendre integration. Polynomials of degree 2*numGaussPoints-1 are integrated
  // exact
  size_t numGaussPoints = (degree + 1) / 2 + numAdditionalPoints;
  base::DataVector roots;
  base::DataVector quadratureweights;
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  double sum = 0.0;

  // on low levels where Lagrange polynomials instead of Bsplines are used the number of Gauss
  // points are not degree dependent and there is only on segment: the whole [0,1] interval
  if ((xValues.size() == 1) || (degree == 3 && (xValues.size() < 5)) ||
      ((degree == 5) && (xValues.size() < 9))) {
    //    std::cout << "Lagrange" << std::endl;
    numGaussPoints = xValues.size() + numAdditionalPoints;
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x = roots[i];
      double basisvalue = expUniformNakBspline(x, degree, index, xValues);
      double integrand = basisvalue * this->weight_function(x);
      sum += integrand * quadratureweights[i];
    }

  } else {
    size_t first_segment = std::max(degree, index);
    size_t last_segment = std::min(xValues.size(), index + degree + 1);
    quadRule.getLevelPointsAndWeightsNormalized(
        std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
    std::vector<double> xi = createNakKnots(xValues, degree);
    for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
      double l = std::max(0.0, xi[segmentIndex]);
      double r = std::min(1.0, xi[segmentIndex + 1]);
      double segmentWidth = r - l;

      for (size_t i = 0; i < roots.getSize(); ++i) {
        // transform roots[i], the quadrature points to the segment on which the Bspline is
        // evaluated and from there to the interval[a,b] on which the weight function is defined
        double x = l + segmentWidth * roots[i];

        double basisvalue = expUniformNakBspline(x, degree, index, xValues);
        double integrand = basisvalue * this->weight_function(x);
        // multiply weights by length_old_interval / length_new_interval
        sum += integrand * quadratureweights[i] * segmentWidth;
      }
    }
  }
  return (b - a) * sum;
}

void BSplineQuadratureEvaluator::calculate1DBSplineIntegrals(
    std::vector<double>& points, std::vector<FloatScalarVector>& basisValues,
    size_t incrementQuadraturePoints, double tol) {
  basisValues.resize(points.size());
  std::vector<FloatScalarVector> newBasisValues(points.size());

  //  std::cout << "BsplineQuadrature calculate1dBsplineIntegrals is called" << std::endl;

  // iteratively increases the numAdditionalPoints until the product of B spline and weight
  // function is exactly inctegrated
  // the numAdditionalPoints of the last index is used as an initial guess for the
  // numAdditionalPoints of the next index. This is serial and must be changed for parallelization
  size_t lastNumAdditionalPoints = 0;
  for (size_t index = 0; index < points.size(); ++index) {
    double err = 1e14;
    // ToDo(rehmemk) optimize this!
    numAdditionalPoints = lastNumAdditionalPoints;
    basisValues[index] = FloatScalarVector(get1DIntegral(points, index));
    if (isCustomWeightFunction) {
      while (err > tol) {
        lastNumAdditionalPoints = numAdditionalPoints;
        numAdditionalPoints += incrementQuadraturePoints;
        // recalculate and check for error < tol
        newBasisValues[index] = FloatScalarVector(get1DIntegral(points, index));
        err = std::fabs(newBasisValues[index].getValue() - basisValues[index].getValue());
        basisValues[index] = newBasisValues[index];

        if (numAdditionalPoints > 490) {
          break;
        }
      }
    }
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
      degree(3),
      a(0),
      b(1) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = 3;
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      numAdditionalPoints(0),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      degree(degree),
      a(0),
      b(1) {
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
      degree(degree),
      a(0),
      b(1) {
  evalConfig.type = CombiEvaluatorTypes::Scalar_BSplineQuadrature;
  evalConfig.degree = degree;
}

BSplineQuadratureEvaluator::BSplineQuadratureEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
    double a, double b, bool normalizeWeights)
    : weight_function(weight_function),
      numAdditionalPoints(numAdditionalPoints),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      degree(degree),
      a(a),
      b(b) {
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
      degree(other.degree),
      a(other.a),
      b(other.b) {
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
