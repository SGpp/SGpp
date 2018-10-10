// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrix.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

namespace sgpp {
namespace optimization {

class ASMatrixGradientMC : public ASMatrix {
 public:
  ASMatrixGradientMC(WrapperScalarFunction objectiveFunc,
                     WrapperScalarFunctionGradient objectiveFuncGradient)
      : ASMatrix(objectiveFunc), objectiveFuncGradient(objectiveFuncGradient) {}
  void createMatrix(size_t numPoints) { createMatrixMonteCarlo(numPoints); };
  void createMatrixMonteCarlo(size_t numPoints);

 private:
  // gradients objective function
  WrapperScalarFunctionGradient objectiveFuncGradient;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
