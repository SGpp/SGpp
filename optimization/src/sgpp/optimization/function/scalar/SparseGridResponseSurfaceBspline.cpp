// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/function/scalar/SparseGridResponseSurfaceBspline.hpp>

#include <algorithm>

namespace sgpp {
namespace optimization {

void SparseGridResponseSurfaceBspline::regular(size_t level) {
  grid->getGenerator().regular(level);
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void SparseGridResponseSurfaceBspline::regularByPoints(size_t numPoints, bool verbose) {
  // set initial level
  size_t level = 1;
  if (boundary == true) level = 0;
  do {
    // todo (rehmemk) instead of trying out until pointnumber matches, use formula for number of
    // grid points
    grid->getStorage().clear();
    grid->getGenerator().regular(level);
    if (verbose == true)
      std::cout << "level " << level << " with " << grid->getSize() << " points\n";
    level++;
  } while (grid->getSize() < numPoints);
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void SparseGridResponseSurfaceBspline::surplusAdaptive(size_t maxNumGridPoints, size_t initialLevel,
                                                       size_t refinementsNum, bool verbose) {
  regular(initialLevel);
  while (grid->getSize() < maxNumGridPoints) {
    refineSurplusAdaptive(refinementsNum);
    interpolant =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
    if (verbose)
      std::cout << "Refining. Calculated a grid with " << grid->getSize() << " points.\n";
  }
}

void SparseGridResponseSurfaceBspline::ritterNovak(size_t maxNumGridPoints, double gamma,
                                                   bool verbose) {
  // this uses the default values for Ritter Novaks initial Level, max level, powerMethod
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(*objectiveFunc, *grid,
                                                                maxNumGridPoints, gamma);
  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
  }
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
  if (verbose) std::cout << "Refining. Calculated a grid with " << grid->getSize() << " points.\n";
}

double SparseGridResponseSurfaceBspline::eval(sgpp::base::DataVector v) {
  transformPoint(v, lb, ub, unitLBounds, unitUBounds);

  return interpolant->eval(v);
}
double SparseGridResponseSurfaceBspline::evalGradient(sgpp::base::DataVector v,
                                                      sgpp::base::DataVector& gradient) {
  transformPoint(v, lb, ub, unitLBounds, unitUBounds);
  double evalValue = interpolantGradient->eval(v, gradient);
  // scale gradient components according to the inner derivative of the chain rule when transforming
  // the interpolation point from the original coordinates to unit cube
  for (size_t i = 0; i < gradient.getSize(); i++) {
    gradient[i] /= (ub[i] - lb[i]);
  }
  return evalValue;
}

double SparseGridResponseSurfaceBspline::getIntegral() {
  sgpp::base::OperationQuadrature* opQuad = sgpp::op_factory::createOperationQuadrature(*grid);
  double unitIntegral = opQuad->doQuadrature(coefficients);
  return unitIntegral * domainVolume();
}

double SparseGridResponseSurfaceBspline::getMean(sgpp::base::DistributionsVector pdfs,
                                                 size_t quadOrder) {
  sgpp::base::OperationWeightedQuadrature* opWQuad =
      sgpp::op_factory::createOperationWeightedQuadrature(*grid);
  mean = opWQuad->doWeightedQuadrature(coefficients, pdfs, quadOrder);
  computedMeanFlag = true;
  return mean;
}

sgpp::base::DataVector SparseGridResponseSurfaceBspline::getVariance(
    sgpp::base::DistributionsVector pdfs, size_t quadOrder) {
  if (!computedMeanFlag) {
    double dummy = getMean(pdfs, quadOrder);
  }
  sgpp::base::OperationWeightedSecondMoment* opWSM =
      sgpp::op_factory::createOperationWeightedSecondMoment(*grid);
  double meanSquare = opWSM->doWeightedQuadrature(coefficients, pdfs, quadOrder);
  std::cout << std::setprecision(16) << "Variance Calculation. Mean: " << mean
            << " meanSquare: " << meanSquare << "\n";
  double variance = meanSquare - mean * mean;
  sgpp::base::DataVector returnVec(3);
  returnVec[0] = variance;
  returnVec[1] = meanSquare;
  returnVec[2] = mean;
  return returnVec;
}

sgpp::base::DataVector SparseGridResponseSurfaceBspline::optimize() {
  /**
   * The gradient method needs a starting point.
   *  We use the grid point with the smallest (most promising) function value and save it in x0.
   */
  sgpp::base::DataVector x0(numDim);
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  // index of grid point with minimal function value
  size_t x0Index =
      std::distance(functionValues.getPointer(),
                    std::min_element(functionValues.getPointer(),
                                     functionValues.getPointer() + functionValues.getSize()));

  x0 = gridStorage.getCoordinates(gridStorage[x0Index]);

  sgpp::optimization::optimizer::GradientDescent gradientDescent(*interpolant,
                                                                 *interpolantGradient);
  gradientDescent.setStartingPoint(x0);
  gradientDescent.optimize();
  const sgpp::base::DataVector& unitXOpt = gradientDescent.getOptimalPoint();
  sgpp::base::DataVector xOpt(unitXOpt);
  transformPoint(xOpt, unitLBounds, unitUBounds, lb, ub);
  return xOpt;
}

// ----------------- auxiliary routines -----------

void SparseGridResponseSurfaceBspline::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  calculateInterpolationCoefficients();
}

void SparseGridResponseSurfaceBspline::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  functionValues.resizeZero(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      p[d] = gridStorage.getPointCoordinate(i, d);
    }
    transformPoint(p, unitLBounds, unitUBounds, lb, ub);
    functionValues[i] = objectiveFunc->eval(p);
  }
  sgpp::base::sle_solver::Armadillo sleSolver;
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "Solving failed!" << std::endl;
  }
}

}  // namespace optimization
}  // namespace sgpp
