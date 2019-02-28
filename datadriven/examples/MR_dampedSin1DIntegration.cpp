// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <iostream>

double func(sgpp::base::DataVector y) {
  // transform to [l,r] = [0, 2.82842708536
  y[0] = 2.82842708536 * y[0];
  return sin(y[0] * sqrt(8) * 0.75 + 1) / (y[0] * sqrt(8) * 0.75 + 1);
}

double l2Error(std::shared_ptr<sgpp::base::NakBsplineBoundaryGrid> grid,
               sgpp::base::DataVector coefficients, sgpp::optimization::WrapperScalarFunction func,
               size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunction I(*grid, coefficients);
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(func.getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    // random points
    // sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector,
    // 0.0, 1.0);
    // uniform 1D points
    randomVector[0] = static_cast<double>(i) / static_cast<double>(numMCPoints);
    double evalInterpolant = I.eval(randomVector);
    double evalObjectiveFunc = func.eval(randomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);

  size_t degree = 3;
  size_t numDim = 1;
  size_t initialLevel = 2;
  size_t numRefine = 50;
  size_t maxNumGridPoints = 2000;
  size_t numMCPoints = 10;
  sgpp::optimization::WrapperScalarFunction Func(numDim, func);
  auto grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
  grid->getGenerator().regular(initialLevel);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::DataVector coefficients;
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gridPointVector[0] = gp.getStandardCoordinate(0);
    functionValues[i] = Func.eval(gridPointVector);
  }
  coefficients.resizeZero(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
  }

  double l2Err = l2Error(grid, coefficients, Func, numMCPoints);
  std::cout << grid->getSize() << " " << l2Err << "\n";

  while (grid->getSize() < maxNumGridPoints) {
    sgpp::base::SurplusRefinementFunctor functor(coefficients, numRefine);
    grid->getGenerator().refine(functor);
    functionValues.resizeZero(gridStorage.getSize());
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
      sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
      gp.getStandardCoordinates(gridPointVector);
      functionValues[i] = Func.eval(gridPointVector);
    }

    coefficients.resizeZero(functionValues.getSize());
    sgpp::optimization::HierarchisationSLE hierSLE(*grid);
    sgpp::optimization::sle_solver::Armadillo sleSolver;
    if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
      std::cout << "ASMatrixNakBspline: Solving failed.\n";
    }
    double l2Err = l2Error(grid, coefficients, Func, numMCPoints);
    std::cout << grid->getSize() << " " << l2Err << "\n";
  }

  return 0;
}
