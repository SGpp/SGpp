// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASResponseSurface.hpp>

namespace sgpp {
namespace optimization {

double ASResponseSurface::l2Error(sgpp::optimization::WrapperScalarFunction objectiveFunc,
                                  size_t numMCPoints) {
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(objectiveFunc.getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    double evalInterpolant = interpolant->eval(randomVector);
    double evalObjectiveFunc = objectiveFunc.eval(randomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
