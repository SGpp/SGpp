// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialScalarProductEvaluator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator()
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      numAdditionalPoints(0) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator(
    sgpp::combigrid::SingleFunction weight_function, bool normalizeWeights,
    size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      numAdditionalPoints(numAdditionalPoints) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::PolynomialScalarProductEvaluator(
    PolynomialScalarProductEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      numAdditionalPoints(other.numAdditionalPoints) {
  evalConfig.type = CombiEvaluatorTypes::Multi_PolynomialScalarProduct;
}

PolynomialScalarProductEvaluator::~PolynomialScalarProductEvaluator() {}

double PolynomialScalarProductEvaluator::quad(LagrangePolynom& p_i, LagrangePolynom& p_j) {
  size_t degree_i = p_i.points.size() - 1;
  size_t degree_j = p_j.points.size() - 1;
  size_t numGaussPoints = (degree_i + degree_j + 3) / 2;

  // performing Gauss-Legendre integration
  auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
  base::DataVector roots;
  base::DataVector quadratureweights;

  // do iterative 1d quadrature
  double scalarProduct_ij = 0.0, scalarProduct_ij_old = 0.0;
  double err = 1e14;
  size_t iteration = 0;
  while (err > 1e-13 && numGaussPoints < quadRule.getMaxSupportedLevel()) {
    quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, quadratureweights);

    scalarProduct_ij = 0.0;
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x_unit = roots[i], w = quadratureweights[i];
      scalarProduct_ij += w * p_i.evaluate(x_unit) * p_j.evaluate(x_unit) * weight_function(x_unit);
    }

    if (iteration > 0) {
      err = std::fabs(scalarProduct_ij_old - scalarProduct_ij);
    }
    scalarProduct_ij_old = scalarProduct_ij;
    numGaussPoints += 1 + numAdditionalPoints;
    iteration += 1;
  }

  return scalarProduct_ij;
}

FloatArrayVector PolynomialScalarProductEvaluator::get1DMixedIntegral(std::vector<double>& points,
                                                                      size_t index_i) {
  FloatArrayVector scalarProducts;

  LagrangePolynom p_i;
  p_i.points = points;
  p_i.point = index_i;

  LagrangePolynom p_j;
  p_j.points = points;

  for (size_t index_j = 0; index_j < points.size(); index_j++) {
    p_j.point = index_j;
    double value = 0.0;

    size_t key = generateKey(index_i, index_j);
    auto it_value = scalarProductsMap.find(key);
    if (it_value != scalarProductsMap.end()) {
      value = it_value->second;
    } else {
      value = quad(p_i, p_j);
      scalarProductsMap[key] = value;
    }

    scalarProducts.at(index_j) = value;
  }

  return scalarProducts;
}

void PolynomialScalarProductEvaluator::calculate1DPolynomialScalarProducts(
    std::vector<double>& points, std::vector<FloatArrayVector>& integrals) {
  for (size_t index_i = 0; index_i < points.size(); ++index_i) {
    integrals.push_back(get1DMixedIntegral(points, index_i));
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

std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector> >
PolynomialScalarProductEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector> >(
      new PolynomialScalarProductEvaluator(*this));
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
