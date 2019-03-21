// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <functional>
#include <iostream>

#include "../../optimization/src/sgpp/optimization/function/scalar/SparseGridResponseSurfaceBspline.hpp"

double f(sgpp::base::DataVector v) { return exp(-v.sum()); }

int main() {
  size_t numDim = 4;
  auto objFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(
      sgpp::optimization::WrapperScalarFunction(numDim, f));

  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineExtended;
  size_t degree = 3;
  size_t numPoints = 5000;
  size_t initialLevel = 2;
  size_t numRefine = 250;

  auto reSurf = sgpp::optimization::SparseGridResponseSurfaceBspline(objFunc, gridType, degree);
  reSurf.surplusAdaptive(numPoints, initialLevel, numRefine);
  std::cout << reSurf.l2Error(objFunc, 10000) << "\n";
}
