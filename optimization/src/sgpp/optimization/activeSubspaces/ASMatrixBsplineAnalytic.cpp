// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include "ASMatrixBsplineAnalytic.hpp"

namespace sgpp {
namespace optimization {

void ASMatrixBsplineAnalytic::buildRegularInterpolant(size_t level) {
  grid->getGenerator().regular(level);
  this->calculateInterpolationCoefficients();
}

void ASMatrixBsplineAnalytic::buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                                       size_t refinementsNum) {
  grid->getGenerator().regular(initialLevel);
  this->calculateInterpolationCoefficients();
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum);
  }
}

void ASMatrixBsplineAnalytic::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  evaluationPoints.resizeZero(gridStorage.getSize(), numDim);
  functionValues.resizeZero(gridStorage.getSize());

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::DataVector gridPointVector = gridStorage.getPointCoordinates(i);
    evaluationPoints.setRow(i, gridPointVector);
    functionValues[i] = objectiveFunc->eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  coefficients = alpha;
}

double ASMatrixBsplineAnalytic::l2InterpolationError(size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunction interpolant(*grid, coefficients);
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(objectiveFunc->getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    double evalInterpolant = interpolant.eval(randomVector);
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

sgpp::base::DataVector ASMatrixBsplineAnalytic::l2InterpolationGradientError(
    std::shared_ptr<sgpp::optimization::WrapperScalarFunctionGradient> objectiveFuncGradient,
    size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunctionGradient interpolant(*grid, coefficients);
  size_t numDim = objectiveFuncGradient->getNumberOfParameters();
  sgpp::base::DataVector errors(numDim, 0);
  sgpp::base::DataVector randomVector(numDim);
  sgpp::base::DataVector interpolantEval(numDim);
  sgpp::base::DataVector gradientEval(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    interpolant.eval(randomVector, interpolantEval);
    objectiveFuncGradient->eval(randomVector, gradientEval);
    for (size_t d = 0; d < numDim; d++) {
      errors[d] += std::pow(interpolantEval[d] - gradientEval[d], 2.0);
    }
  }
  for (size_t d = 0; d < numDim; d++) {
    errors[d] = sqrt(errors[d] / static_cast<double>(numMCPoints));
  }
  return errors;
}

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
