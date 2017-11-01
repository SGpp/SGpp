// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator_1d.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

BSplineInterpolationEvaluator_1d::BSplineInterpolationEvaluator_1d()
    : evaluationPoint(0.0), basisValues(), xValues(), degree(3) {}

BSplineInterpolationEvaluator_1d::BSplineInterpolationEvaluator_1d(size_t degree)
    : evaluationPoint(0.0), basisValues(), xValues(), degree(degree) {}

BSplineInterpolationEvaluator_1d::~BSplineInterpolationEvaluator_1d() {}

BSplineInterpolationEvaluator_1d::BSplineInterpolationEvaluator_1d(
    const BSplineInterpolationEvaluator_1d& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues),
      degree(other.degree) {}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
BSplineInterpolationEvaluator_1d::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>(
      new BSplineInterpolationEvaluator_1d(*this));
}

bool BSplineInterpolationEvaluator_1d::needsOrderedPoints() { return true; }

bool BSplineInterpolationEvaluator_1d::needsParameter() { return true; }

void BSplineInterpolationEvaluator_1d::setDegree(size_t const& deg) { degree = deg; }

void BSplineInterpolationEvaluator_1d::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisValues();
}

void BSplineInterpolationEvaluator_1d::setGridPoints(std::vector<double> const& x) {
  xValues = x;
  computeBasisValues();
}

void BSplineInterpolationEvaluator_1d::computeBasisValues() {
  /*
   * Am Rand: Spiegle Gitterweiten nach außen
   *
   * nonuniform B-Spline vom Grad n:
   * Knotenfolge xi_k,.,xi_{k+n+1}
   * b^n_{k,xi} = gamma b^{n-1}_{k,xi} + (1-gamma) b^{n-1}_{k+1,xi}
   * gamma(x) = (x-xi_k)/(xi_{k+n}-xi_k)
   *
   *
   */

  basisValues.resize(xValues.size(), sgpp::combigrid::FloatScalarVector(0));

  if (xValues.size() == 1) {
    basisValues[0] = 1.0;
    return;
  }
  // Lagrange polynomials for less than 9 points because 9 is the number of gridpoints of a uniform
  // boundary grid of level 3 and this is the first level with enough gridpoints for nak B-Splines
  // Should work for degree 5 as well
  // For degree 7 and higher level 3 with nak is too small to provide enough knots even for one
  // single spline
  else if (xValues.size() < 9) {
    for (size_t i = 0; i < xValues.size(); i++) {
      basisValues[i] = LagrangePolynomial(evaluationPoint, xValues, i);
    }
    return;
  }
  std::vector<double> xi;
  createKnots(xValues, degree, xi);
  for (size_t i = 0; i < xValues.size(); i++) {
    basisValues[i] = nonUniformBSpline(evaluationPoint, degree, i, xi);
  }
}

void BSplineInterpolationEvaluator_1d::setFunctionValuesAtGridPoints(
    std::vector<double>& pFunctionValues) {
  size_t numGridPoints=xValues.size();

  std::cout<<"test";


  sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
  sgpp::base::DataVector coefficients_sle(numGridPoints);
  sgpp::base::DataVector functionValues(numGridPoints);

  for (size_t ixEvalPoints = 0; ixEvalPoints<xValues.size(); ++ixEvalPoints ) {
    
    functionValues[ixEvalPoints] = pFunctionValues[ixEvalPoints];

    std::vector<double> basisValues_vec(basisValues.size());
    for (size_t i = 0; i < basisValues.size(); i++) {
      basisValues_vec[i] = basisValues[i].value();
    }

   
    // iterate over every gridpoint
    for (size_t ixBasisFunctions = 0;ixBasisFunctions < basisValues_vec.size() ; ++ixBasisFunctions ) {
      double splineValue = 1.0;


      splineValue *= basisValues_vec[ixBasisFunctions];
      A.set(ixEvalPoints, ixBasisFunctions, splineValue);
    }
    
  }

  sgpp::optimization::FullSLE sle(A);
  sgpp::optimization::sle_solver::Auto solver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  bool solved = solver.solve(sle, functionValues, coefficients_sle);

  std::cout << A.toString() << std::endl;

  basisCoefficients = coefficients_sle;
}

} /* namespace combigrid */
} /* namespace sgpp */
