// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {

// ToDo (rehmemk) This is a lot less efficient than LTwoScalarProductHashMapNakBsplineBoundary.
// Probably  because of the missing hash map feature
FloatArrayVector BSplineScalarProductEvaluator::get1DL2ScalarProduct(
    std::vector<double> const& points, size_t const& index_i) {
  FloatArrayVector sums;

  for (size_t index_j = 0; index_j < points.size(); index_j++) {
    // performing Gauss-Legendre integration with twice as many points as for the simple integrals
    size_t numGaussPoints = 2 * ((degree + 1) / 2 + numAdditionalPoints);
    sgpp::base::DataVector roots;
    sgpp::base::DataVector quadratureweights;
    auto& quadRule = base::GaussLegendreQuadRule1D::getInstance();
    double sum = 0.0;

    // constant function for single point, Lagrange polynomials while not enough knots for not a
    // knot B-splines, nak B-splines otherwise
    if ((xValues.size() == 1) || (degree == 3 && (xValues.size() < 5)) ||
        ((degree == 5) && (xValues.size() < 9))) {
      numGaussPoints = 2 * xValues.size() + numAdditionalPoints;
      quadRule.getLevelPointsAndWeightsNormalized(
          std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x = roots[i];
        double productValue = expUniformNakBspline(x, degree, index_i, xValues) *
                              expUniformNakBspline(x, degree, index_j, xValues);
        double integrand = productValue * this->weight_function(x);
        sum += integrand * quadratureweights[i];
      }
    } else {
      // only if supports of B-splines overlap the scalar product must be calculated
      if (std::abs(static_cast<int>(index_i) - static_cast<int>(index_j)) <=
          static_cast<int>(degree)) {
        quadRule.getLevelPointsAndWeightsNormalized(
            std::min(numGaussPoints, quadRule.getMaxSupportedLevel()), roots, quadratureweights);
        size_t first_segment_j = std::max(degree, index_j);
        size_t last_segment_j = std::min(xValues.size(), index_j + degree + 1);
        size_t first_segment_i = std::max(degree, index_i);
        size_t last_segment_i = std::min(xValues.size(), index_i + degree + 1);
        size_t first_segment = std::min(first_segment_i, first_segment_j);
        size_t last_segment = std::max(last_segment_i, last_segment_j);
        // ToDO(rehmemk) Rewrite this whole  routine , don't use createKnots
        std::vector<double> xi = createNakKnots(xValues, degree);
        for (size_t segmentIndex = first_segment; segmentIndex < last_segment; segmentIndex++) {
          double l = std::max(0.0, xi[segmentIndex]);
          double r = std::min(1.0, xi[segmentIndex + 1]);
          double segmentWidth = r - l;

          for (size_t i = 0; i < roots.getSize(); ++i) {
            double x = l + segmentWidth * roots[i];
            double productValue = expUniformNakBspline(x, degree, index_i, xValues) *
                                  expUniformNakBspline(x, degree, index_j, xValues);
            double integrand = productValue * this->weight_function(x);
            // multiply weights by length_old_interval / length_new_interval
            sum += integrand * quadratureweights[i] * segmentWidth;
          }
        }
      }
    }
    sums[index_j] = (b - a) * sum;
  }
  return sums;
}

void BSplineScalarProductEvaluator::calculate1DBSplineScalarProducts(
    std::vector<double>& points, std::vector<FloatArrayVector>& basisValues,
    size_t incrementQuadraturePoints, double tol) {
  basisValues.resize(points.size());
  std::vector<FloatArrayVector> newBasisValues(points.size());

  // iteratively increases the numAdditionalPoints until the product of B spline and weight
  // function is exactly inctegrated
  // the numAdditionalPoints of the last index is used as an initial guess for the
  // numAdditionalPoints of the next index. This is serial and must be changed for parallelization
  size_t lastNumAdditionalPoints = 0;

  // #pragma omp parallel for schedule(dynamic)
  for (size_t index_i = 0; index_i < points.size(); ++index_i) {
    double err = 1e14;
    numAdditionalPoints = lastNumAdditionalPoints;
    basisValues[index_i] = get1DL2ScalarProduct(points, index_i);

    // ToDo (rehmemk) This is always false!!!
    if (isCustomWeightFunction) {
      while (err > tol) {
        lastNumAdditionalPoints = numAdditionalPoints;
        numAdditionalPoints += incrementQuadraturePoints;
        // recalculate and check for error < tol
        newBasisValues[index_i] = get1DL2ScalarProduct(points, index_i);
        basisValues[index_i].sub(newBasisValues[index_i]);
        err = basisValues[index_i].norm();
        basisValues[index_i] = newBasisValues[index_i];

        //        std::cout << numAdditionalPoints << " " << err << std::endl;

        if (numAdditionalPoints > 490) {
          break;
        }
      }
    }
    //    std::cout << "index i: " << index_i << std::endl;
    //    for (size_t j = 0; j < basisValues[index_i].size(); j++) {
    //      std::cout << basisValues[index_i].get(j) << " ";
    //    }
    //    std::cout << "\n";
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
      degree(3),
      a(0),
      b(1),
      numAdditionalPoints(0) {
  evalConfig.type = CombiEvaluatorTypes::Multi_BSplineScalarProduct;
  evalConfig.degree = 3;
}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(size_t degree)
    : weight_function(constantFunction<double>(1.0)),
      normalizeWeights(false),
      isCustomWeightFunction(false),
      degree(degree),
      a(0),
      b(1),
      numAdditionalPoints(0) {
  evalConfig.type = CombiEvaluatorTypes::Multi_BSplineScalarProduct;
  evalConfig.degree = degree;
}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(
    size_t degree, sgpp::combigrid::SingleFunction weight_function, double a, double b,
    bool normalizeWeights, size_t numAdditionalPoints)
    : weight_function(weight_function),
      normalizeWeights(normalizeWeights),
      isCustomWeightFunction(true),
      degree(degree),
      a(a),
      b(b),
      numAdditionalPoints(numAdditionalPoints) {
  evalConfig.type = CombiEvaluatorTypes::Multi_BSplineScalarProduct;
  evalConfig.degree = degree;
}

BSplineScalarProductEvaluator::BSplineScalarProductEvaluator(
    BSplineScalarProductEvaluator const& other)
    : xValues(other.xValues),
      basisValues(other.basisValues),
      weight_function(other.weight_function),
      normalizeWeights(other.normalizeWeights),
      isCustomWeightFunction(other.isCustomWeightFunction),
      degree(other.degree),
      a(other.a),
      b(other.b),
      numAdditionalPoints(other.numAdditionalPoints) {
  evalConfig.type = CombiEvaluatorTypes::Multi_BSplineScalarProduct;
  evalConfig.degree = other.degree;
}

void BSplineScalarProductEvaluator::setParameter(const FloatArrayVector& param) { return; }

void BSplineScalarProductEvaluator::setBasisCoefficientsAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

} /* namespace combigrid */
} /* namespace sgpp*/
