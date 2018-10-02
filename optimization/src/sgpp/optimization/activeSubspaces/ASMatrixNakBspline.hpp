// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrix.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

namespace sgpp {
namespace optimization {

class ASMatrixNakBspline : public ASMatrix {
 public:
  ASMatrixNakBspline(WrapperScalarFunction objectiveFunc, size_t degree)
      : ASMatrix(objectiveFunc), degree(degree) {}
  void buildRegularInterpolant(size_t level);
  void createMatrix(size_t numPoints);
  void createMatrixMonteCarlo(size_t numPoints);

  // ToDo (rehmemk) createMatrix routine using analytical integrals of the B-Splines

 private:
  size_t degree;
  sgpp::base::DataVector coefficients;
  std::shared_ptr<sgpp::base::Grid> grid;
  bool interpolantFlag = 0;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
