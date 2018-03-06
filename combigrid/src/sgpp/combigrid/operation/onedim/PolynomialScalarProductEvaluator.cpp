// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialScalarProductEvaluator.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace combigrid {

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0),
      xlower(0.0),
      xupper(1.0) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator(
    sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints),
      xlower(0.0),
      xupper(1.0) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator(
    PolynomialScalarProductEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints),
      xlower(other.xlower),
      xupper(other.xupper) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> orthogBasis)
    : weight_function([orthogBasis](double x) { return orthogBasis->pdf(x); }),
      normalizeWeights(false),
      isCustomWeightFunction(true),
      numAdditionalPoints(0),
      xlower(0.0),
      xupper(1.0) {
  initializeBounds(orthogBasis);
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::~PolynomialScalarProductEvaluator() {}

void PolynomialScalarProductEvaluator::initializeBounds(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> orthogBasis) {
#ifdef USE_DAKOTA
  Pecos::RealRealPair bounds = orthogBasis->getRandomVariable()->bounds();
  xlower = bounds.first;
  xupper = bounds.second;
#endif
}

double PolynomialScalarProductEvaluator::quad(LagrangePolynom& p_i, LagrangePolynom& p_j) {
  size_t degree_i = p_i.degree() - 1;
  size_t degree_j = p_j.degree() - 1;
  size_t numGaussPoints = std::max(static_cast<size_t>(1), (degree_i + degree_j + 3) / 2);

  // performing Gauss-Legendre integration
  auto func = [&p_i, &p_j, this](double x_unit, double x_prob) {
    return p_i.evaluate(x_unit) * p_j.evaluate(x_unit) * weight_function(x_prob);
  };

  return GaussLegendreQuadrature(numGaussPoints)
      .evaluate_iteratively(func, xlower, xupper, numGaussPoints, 1 + numAdditionalPoints);
}

FloatArrayVector PolynomialScalarProductEvaluator::get1DMixedIntegral(std::vector<double>& points,
                                                                      size_t index_i) {
  FloatArrayVector scalarProducts;

  LagrangePolynom p_i;
  p_i.initialize(index_i, points);

  LagrangePolynom p_j;
  for (size_t index_j = 0; index_j < points.size(); index_j++) {
    p_j.initialize(index_j, points);
    scalarProducts.at(index_j) = quad(p_i, p_j);
  }

  return scalarProducts;
}

void PolynomialScalarProductEvaluator::calculate1DPolynomialScalarProducts(
    std::vector<double>& points, std::vector<FloatArrayVector>& integrals) {
  basisValues.resize(points.size());
  // #pragma omp parallel for schedule(static)
  for (size_t index_i = 0; index_i < points.size(); ++index_i) {
    integrals[index_i] = get1DMixedIntegral(points, index_i);
  }
}

bool PolynomialScalarProductEvaluator::needsOrderedPoints() { return false; }

bool PolynomialScalarProductEvaluator::needsParameter() { return false; }

void PolynomialScalarProductEvaluator::setGridPoints(std::vector<double> const& newXValues) {
  xValues = newXValues;
  basisValues.clear();
  calculate1DPolynomialScalarProducts(xValues, basisValues);

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

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>
PolynomialScalarProductEvaluator::cloneLinear() {
  return std::make_shared<PolynomialScalarProductEvaluator>(*this);
}

void PolynomialScalarProductEvaluator::setParameter(const FloatArrayVector& param) { return; }

void PolynomialScalarProductEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

size_t PolynomialScalarProductEvaluator::generateKey(size_t i, size_t j) {
  // sort for commuativity
  size_t smaller_index = 0;
  size_t larger_index = 0;
  if (i < j) {
    smaller_index = i;
    larger_index = j;
  } else {
    smaller_index = j;
    larger_index = i;
  }

  // cantor pairing function
  return (smaller_index + larger_index) * (smaller_index + larger_index + 1) / 2 + larger_index;
}

} /* namespace combigrid */
} /* namespace sgpp*/
