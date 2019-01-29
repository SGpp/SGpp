// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/datadriven/activeSubspaces/ResponseSurface.hpp>

namespace sgpp {
namespace datadriven {

double ResponseSurface::l2Error(std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
                                size_t numMCPoints) {
  double l2Err = 0.0;
  sgpp::base::DataVector uniformRandomVector(objectiveFunc->getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(uniformRandomVector, 0.0,
                                                                          1.0);
    double evalInterpolant = this->eval(uniformRandomVector);
    double evalObjectiveFunc = objectiveFunc->eval(uniformRandomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
