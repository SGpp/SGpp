// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/optimization/function/scalar/ResponseSurface.hpp>

namespace sgpp {
namespace optimization {

double ResponseSurface::eval(sgpp::base::DataVector v) { return interpolant->eval(v); }

double ResponseSurface::l2Error(std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
                                size_t numMCPoints) {
  double l2Err = 0.0;
  size_t numDim = objectiveFunc->getNumberOfParameters();
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] =
          sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    double evalInterpolant = this->eval(randomVector);
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);

    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));

  return l2Err;
}

sgpp::base::DataVector ResponseSurface::nrmsError(
    std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc, size_t numMCPoints) {
  double max = 0.0;
  double min = 100000.0;
  double l2Err = 0.0;
  size_t numDim = objectiveFunc->getNumberOfParameters();
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] =
          sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
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
  double vol = 0;
  for (size_t d = 0; d < numDim; d++) {
    vol += ub[d] - lb[d];
  }
  return vol;
}

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
