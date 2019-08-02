// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp>

#include <iostream>

double func1(sgpp::base::DataVector v) { return 1; }
double funcx(sgpp::base::DataVector v) { return v[0]; }
double funcx2(sgpp::base::DataVector v) { return v[0] * v[0]; }
double funcx3(sgpp::base::DataVector v) { return v[0] * v[0] * v[0]; }
double funcx4(sgpp::base::DataVector v) { return v[0] * v[0] * v[0] * v[0]; }
double funcx5(sgpp::base::DataVector v) { return v[0] * v[0] * v[0] * v[0] * v[0]; }

double l2Error(std::shared_ptr<sgpp::base::NakBsplineExtendedGrid> grid,
               sgpp::base::DataVector coefficients, sgpp::base::WrapperScalarFunction func,
               size_t numMCPoints) {
  sgpp::base::InterpolantScalarFunction I(*grid, coefficients);
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(func.getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    // random points
    // sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(randomVector,
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

double interpolateAndError(size_t degree, size_t dim, size_t level,
                           sgpp::base::WrapperScalarFunction func, size_t numMCPoints) {
  auto grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(dim, degree);
  grid->getGenerator().regular(level);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    f_values[i] = func.eval(gridPointVector);
  }
  sgpp::base::DataVector alpha(gridStorage.getSize());
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  sgpp::base::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed.\n";
    return 0;
  }

  double l2Err = l2Error(grid, alpha, func, numMCPoints);
  return l2Err;
}

sgpp::base::DataVector polynomialError(size_t degree, size_t dim, size_t level,
                                       size_t numMCPoints) {
  sgpp::base::DataVector errors(6);
  sgpp::base::WrapperScalarFunction Func1(dim, func1);
  errors[0] = interpolateAndError(degree, dim, level, Func1, numMCPoints);
  sgpp::base::WrapperScalarFunction Funcx(dim, funcx);
  errors[1] = interpolateAndError(degree, dim, level, Funcx, numMCPoints);
  sgpp::base::WrapperScalarFunction Funcx2(dim, funcx2);
  errors[2] = interpolateAndError(degree, dim, level, Funcx2, numMCPoints);
  sgpp::base::WrapperScalarFunction Funcx3(dim, funcx3);
  errors[3] = interpolateAndError(degree, dim, level, Funcx3, numMCPoints);
  sgpp::base::WrapperScalarFunction Funcx4(dim, funcx4);
  errors[4] = interpolateAndError(degree, dim, level, Funcx4, numMCPoints);
  sgpp::base::WrapperScalarFunction Funcx5(dim, funcx5);
  errors[5] = interpolateAndError(degree, dim, level, Funcx5, numMCPoints);

  return errors;
}

double objectiveFunc(sgpp::base::DataVector v) { return sin(10 * v[0]); }
int main() {
  //  sgpp::base::Printer::getInstance().setVerbosity(-1);
  //  size_t dim = 1;
  //  size_t degree = 5;
  //  size_t numMCPoints = 1000;
  //  for (size_t level = 1; level < 6; level++) {
  //    sgpp::base::DataVector errors = polynomialError(degree, dim, level, numMCPoints);
  //    std::cout << level << ": " << errors.toString() << "\n";
  //  }

  size_t degree = 5;
  size_t numDim = 1;
  auto grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
  grid->getGenerator().regular(5);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    functionValues[i] = objectiveFunc(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector coefficients(functionValues.getSize());
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  sgpp::base::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
  }

  sgpp::base::WrapperScalarFunction func(numDim, objectiveFunc);
  size_t numMCPoints = 10000;
  double l2Err = l2Error(grid, coefficients, func, numMCPoints);

  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineExtended;
  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0)) * 2;
  sgpp::datadriven::NakBsplineScalarProducts scalarProducts(gridType, gridType, degree, degree,
                                                            quadOrder);

  // trick: the first basis function is the constant one function. So the scalar product with this
  // is the integral!
  sgpp::base::DataVector quadCoefficients(coefficients.getSize(), 0.0);
  quadCoefficients[0] = 1;
  double integral =
      scalarProducts.calculateScalarProduct(grid, coefficients, grid, quadCoefficients);
  std::cout << "l2 Error: " << l2Err << "\n";
  std::cout << "integral: " << integral << "\n";

  sgpp::base::SNakBsplineExtendedBase nakBBase(degree);
  double integral2 = 0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    integral2 += coefficients[i] * nakBBase.getIntegral(gridStorage.getPointLevel(i, 0),
                                                        gridStorage.getPointIndex(i, 0));
  }
  std::cout << "integral2: " << integral2 << "\n";

  return 0;
}
