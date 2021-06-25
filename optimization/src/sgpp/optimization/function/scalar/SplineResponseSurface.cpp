// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <sgpp/optimization/function/scalar/SplineResponseSurface.hpp>

namespace sgpp {
namespace optimization {

void SplineResponseSurface::regular(size_t level) {
  grid->getGenerator().regular(level);
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);

  // ToDo (rehmemk) poly grid types don't have gradients. They should go into a seperate class
  if ((gridType != sgpp::base::GridType::PolyBoundary) &&
      (gridType != sgpp::base::GridType::ModPoly)) {
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }
}

void SplineResponseSurface::regularByPoints(size_t numPoints, bool verbose) {
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
  // ToDo (rehmemk) poly grid types don't have gradients. They should go into a seperate class
  if ((gridType != sgpp::base::GridType::PolyBoundary) &&
      (gridType != sgpp::base::GridType::ModPoly)) {
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }
}

void SplineResponseSurface::surplusAdaptive(size_t maxNumGridPoints, size_t initialLevel,
                                            size_t refinementsNum, bool verbose) {
  regular(initialLevel);
  while (grid->getSize() < maxNumGridPoints) {
    refineSurplusAdaptive(refinementsNum, verbose);
    if (verbose)
      std::cout << "Refining. Calculated a grid with " << grid->getSize() << " points.\n";
  }
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  // ToDo (rehmemk) poly grid types don't have gradients. They should go into a seperate class
  if ((gridType != sgpp::base::GridType::PolyBoundary) &&
      (gridType != sgpp::base::GridType::ModPoly)) {
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }
}

void SplineResponseSurface::refineSurplusAdaptive(size_t refinementsNum, bool verbose) {
  if (computedCoefficientsFlag == true) {
    nextSurplusAdaptiveGrid(refinementsNum, verbose);
    if (verbose) std::cout << "Refined to a grid with " << grid->getSize() << " points.\n";
  }
  calculateInterpolationCoefficients();
}

/**
 * refines the grid surplus adaptive but does not recalculate interpolation coefficients
 *@param refinementsNum	number of grid points which should be refined
 *@param verbose        print information on the refine points
 */
void SplineResponseSurface::nextSurplusAdaptiveGrid(size_t refinementsNum, bool verbose) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  computedCoefficientsFlag = false;
}

void SplineResponseSurface::ritterNovak(size_t maxNumGridPoints, double gamma, size_t initialLevel,
                                        bool verbose) {
  // this uses the default values for Ritter Novaks initial Level, max level, powerMethod
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(
      *objectiveFunc, *grid, maxNumGridPoints, gamma, static_cast<base::level_t>(initialLevel));
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

double SplineResponseSurface::eval(sgpp::base::DataVector v) {
  transformPoint(v, lb, ub, unitLBounds, unitUBounds);

  return interpolant->eval(v);
}
double SplineResponseSurface::evalGradient(sgpp::base::DataVector v,
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

double SplineResponseSurface::getIntegral() {
  sgpp::base::OperationQuadrature* opQuad = sgpp::op_factory::createOperationQuadrature(*grid);
  double unitIntegral = opQuad->doQuadrature(coefficients);
  // std::cout << "domainVolume " << domainVolume() << "\n";
  return unitIntegral * domainVolume();
}

double SplineResponseSurface::getMean(sgpp::base::DistributionsVector pdfs, size_t quadOrder) {
  sgpp::base::OperationWeightedQuadrature* opWQuad =
      sgpp::op_factory::createOperationWeightedQuadrature(*grid, quadOrder);
  mean = opWQuad->doWeightedQuadrature(coefficients, pdfs);
  computedMeanFlag = true;
  return mean;
}

sgpp::base::DataVector SplineResponseSurface::getVariance(sgpp::base::DistributionsVector pdfs,
                                                          size_t quadOrder) {
  if (!computedMeanFlag) {
    getMean(pdfs, quadOrder);
  }
  sgpp::base::OperationWeightedSecondMoment* opWSM =
      sgpp::op_factory::createOperationWeightedSecondMoment(*grid, quadOrder);
  double meanSquare = opWSM->doWeightedQuadrature(coefficients, pdfs);
  // std::cout << std::setprecision(16) << "Variance Calculation. Mean: " << mean
  //           << " meanSquare: " << meanSquare << "\n";
  double variance = meanSquare - mean * mean;
  sgpp::base::DataVector returnVec(3);
  returnVec[0] = variance;
  returnVec[1] = meanSquare;
  returnVec[2] = mean;
  return returnVec;
}

sgpp::base::DataVector SplineResponseSurface::optimize() {
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

  sgpp::optimization::optimizer::GradientDescent optimizer(*interpolant, *interpolantGradient);
  // sgpp::optimization::optimizer::AdaptiveGradientDescent optimizer(*interpolant,
  //                                                                  *interpolantGradient);

  optimizer.setStartingPoint(x0);
  optimizer.optimize();
  const sgpp::base::DataVector& unitXOpt = optimizer.getOptimalPoint();
  sgpp::base::DataVector xOpt(unitXOpt);
  transformPoint(xOpt, unitLBounds, unitUBounds, lb, ub);
  return xOpt;
}

// ----------------- auxiliary routines -----------

void SplineResponseSurface::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  functionValues.resizeZero(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    // sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    // for (size_t d = 0; d < gridStorage.getDimension(); d++) {
    //   p[d] = gridStorage.getPointCoordinate(i, d);
    // }
    sgpp::base::DataVector p = gridStorage.getPointCoordinates(i);
    transformPoint(p, unitLBounds, unitUBounds, lb, ub);
    functionValues[i] = objectiveFunc->eval(p);
  }

  if ((gridType == sgpp::base::GridType::PolyBoundary) ||
      (gridType == sgpp::base::GridType::ModPoly)) {
    /**
     * For Bungartz Polynomials use hierarchization operation
     **/
    std::unique_ptr<sgpp::base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(functionValues);
    coefficients = functionValues;
  } else {
    /**
     * For spline basis use SLE
     **/
    std::unique_ptr<sgpp::base::sle_solver::SLESolver> sleSolver;
#ifdef USE_EIGEN
    // sgpp::base::sle_solver::Eigen sleSolver;
    sleSolver.reset(new sgpp::base::sle_solver::Eigen());
#else
    // if compiled without Eigen use simple Gaussian Elimination
    // sgpp::base::sle_solver::GaussianElimination sleSolver;
    sleSolver.reset(new sgpp::base::sle_solver::GaussianElimination());
#endif /* USE_EIGEN */
    sgpp::base::Printer::getInstance().setVerbosity(-1);
    sgpp::base::HierarchisationSLE hierSLE(*grid);
    if (!sleSolver->solve(hierSLE, functionValues, coefficients)) {
      std::cout << "Solving failed!" << std::endl;
    }
  }
  computedCoefficientsFlag = true;
}

}  // namespace optimization
}  // namespace sgpp
