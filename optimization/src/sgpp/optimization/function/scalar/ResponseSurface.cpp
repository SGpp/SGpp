// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/function/scalar/ResponseSurface.hpp>

#include <string>

namespace sgpp {
namespace optimization {

double ResponseSurface::eval(sgpp::base::DataVector v) { return interpolant->eval(v); }

double ResponseSurface::l2Error(std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc,
                                size_t numMCPoints) {
  double l2Err = 0.0;
  size_t numDim = objectiveFunc->getNumberOfParameters();
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = sgpp::base::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    double evalInterpolant = this->eval(randomVector);
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);

    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));

  return l2Err;
}

sgpp::base::DataVector ResponseSurface::nrmsError(
    std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc, size_t numMCPoints) {
  double max = -1e+15;
  double min = 1e+15;
  double l2Err = 0.0;
  size_t numDim = objectiveFunc->getNumberOfParameters();
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = sgpp::base::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    double evalInterpolant = this->eval(randomVector);
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);

    if (evalObjectiveFunc > max) {
      max = evalObjectiveFunc;
    }

    if (evalObjectiveFunc < min) {
      min = evalObjectiveFunc;
    }

    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));

  sgpp::base::DataVector errorVector(4);
  errorVector[0] = l2Err / (max - min);
  errorVector[1] = l2Err;
  errorVector[2] = min;
  errorVector[3] = max;
  return errorVector;
}

sgpp::base::DataVector ResponseSurface::nrmsErrorFromTestData(const std::string& fileName,
                                                              size_t numMCPoints, size_t numDim) {
  sgpp::base::DataMatrix data = sgpp::base::DataMatrix::fromFile(fileName);
  double max = -1e+15;
  double min = 1e+15;
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = data.get(i, d);
    }
    double evalInterpolant = this->eval(randomVector);
    double evalObjectiveFunc = data.get(i, numDim);

    if (evalObjectiveFunc > max) {
      max = evalObjectiveFunc;
    }

    if (evalObjectiveFunc < min) {
      min = evalObjectiveFunc;
    }

    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));

  sgpp::base::DataVector errorVector(4);
  errorVector[0] = l2Err / (max - min);
  errorVector[1] = l2Err;
  errorVector[2] = min;
  errorVector[3] = max;
  return errorVector;
}

void ResponseSurface::precalculateErrorTestData(
    std::shared_ptr<sgpp::base::ScalarFunction> objectiveFunc, size_t numMCPoints,
    const std::string& fileName) {
  size_t numDim = objectiveFunc->getNumberOfParameters();
  sgpp::base::DataMatrix data(numMCPoints, numDim + 1);
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = sgpp::base::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);
    for (size_t d = 0; d < numDim; d++) {
      data.set(i, d, randomVector[d]);
      data.set(i, numDim, evalObjectiveFunc);
    }
  }
  data.toFile(fileName);
}
void ResponseSurface::transformPoint(sgpp::base::DataVector& v, sgpp::base::DataVector lBounds,
                                     sgpp::base::DataVector uBounds,
                                     sgpp::base::DataVector newlBounds,
                                     sgpp::base::DataVector newuBounds) {
  v.sub(lBounds);
  uBounds.sub(lBounds);
  v.componentwise_div(uBounds);
  newuBounds.sub(newlBounds);
  v.componentwise_mult(newuBounds);
  v.add(newlBounds);
}

double ResponseSurface::domainVolume() {
  double vol = 1;
  for (size_t d = 0; d < numDim; d++) {
    vol *= ub[d] - lb[d];
  }
  return vol;
}

}  // namespace optimization
}  // namespace sgpp
