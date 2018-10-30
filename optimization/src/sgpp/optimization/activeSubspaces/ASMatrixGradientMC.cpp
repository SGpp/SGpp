// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>

namespace sgpp {
namespace optimization {

void ASMatrixGradientMC::createMatrixMonteCarlo(
    size_t numPoints, WrapperScalarFunctionGradient objectiveFuncGradient) {
  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  evaluationPoints.resizeZero(numPoints, numDim);
  functionValues.resizeZero(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    sgpp::base::DataVector gradient(numDim);
    functionValues[i] = objectiveFuncGradient.eval(randomVector, gradient);
    evaluationPoints.setRow(i, randomVector);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

void ASMatrixGradientMC::createMatrixMonteCarloFiniteDifference(size_t numPoints, double h) {
  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  evaluationPoints.resizeZero(numPoints, numDim);
  functionValues.resizeZero(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    sgpp::base::DataVector gradient(numDim);
    functionValues[i] = objectiveFunc->eval(randomVector);
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] += h;
      gradient[d] = (objectiveFunc->eval(randomVector) - functionValues[i]) / h;
    }
    evaluationPoints.setRow(i, randomVector);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
