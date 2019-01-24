// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/datadriven/activeSubspaces/ASMatrix.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Used to create, store and use the matrix C for the detection of active subspaces using a Monte
 * Carlo approximation. This is more or less obsolete because we use Paul G Constantines python
 * framework for comparison
 */
class ASMatrixGradientMC : public ASMatrix {
 public:
  explicit ASMatrixGradientMC(std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc)
      : ASMatrix(), objectiveFunc(objectiveFunc), numDim(objectiveFunc->getNumberOfParameters()) {}
  void createMatrix(size_t numPoints) { createMatrixMonteCarloFiniteDifference(numPoints); }
  void createMatrixMonteCarlo(
      size_t numPoints, sgpp::optimization::WrapperScalarFunctionGradient objectiveFuncGradient);
  void createMatrixMonteCarloFiniteDifference(size_t numPoints, double h = 1e-08);

 private:
  // objective function
  std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc;
  // dimensionality
  size_t numDim;
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
