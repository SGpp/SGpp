// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/function/scalar/SparseGridResponseSurfaceBspline.hpp>

#include <algorithm>

namespace sgpp {
namespace optimization {

void SparseGridResponseSurfaceBspline::initialize() {
  numDim = objectiveFunc->getNumberOfParameters();
  if (gridType == sgpp::base::GridType::Bspline) {
    grid = std::make_shared<sgpp::base::BsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SBsplineBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::BsplineBoundary) {
    grid = std::make_shared<sgpp::base::BsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SBsplineBoundaryBase>(degree);
    boundary = true;
  } else if (gridType == sgpp::base::GridType::ModBspline) {
    grid = std::make_shared<sgpp::base::ModBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SBsplineModifiedBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::BsplineClenshawCurtis) {
    grid = std::make_shared<sgpp::base::BsplineClenshawCurtisGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SBsplineClenshawCurtisBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::FundamentalSpline) {
    grid = std::make_shared<sgpp::base::FundamentalSplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SFundamentalSplineBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::ModFundamentalSpline) {
    grid = std::make_shared<sgpp::base::ModFundamentalSplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SFundamentalSplineModifiedBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
    boundary = true;
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
    boundary = false;
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
    boundary = false;
  } else {
    throw sgpp::base::generation_exception(
        "SparseGridResponseSurfaceBspline: gridType not supported.");
  }
}

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

// void SparseGridResponseSurfaceBspline::regularData(size_t level,
//                                                   sgpp::base::DataMatrix evaluationPoints,
//                                                   sgpp::base::DataVector functionValues,
//                                                   double lambda) {
//  grid->getGenerator().regular(level);
//  double mse = 0;
//  sgpp::base::DataVector errorPerBasis;
//  coefficients = sgpp::datadriven::EigenRegression(
//      grid, degree, sgpp::datadrivenDataMatrixToEigen(evaluationPoints), functionValues, mse,
//      errorPerBasis);
//  interpolant =
//      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
//  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
//      *grid, coefficients);
//}

double SparseGridResponseSurfaceBspline::eval(sgpp::base::DataVector v) {
  return interpolant->eval(v);
}
double SparseGridResponseSurfaceBspline::evalGradient(sgpp::base::DataVector v,
                                                      sgpp::base::DataVector& gradient) {
  return interpolantGradient->eval(v, gradient);
}

double SparseGridResponseSurfaceBspline::evalNonUniform(sgpp::base::DataVector v,
                                                        sgpp::base::DataVector lBounds,
                                                        sgpp::base::DataVector uBounds) {
  sgpp::base::DataVector newlBounds(lBounds.getSize(), 0.0);
  sgpp::base::DataVector newuBounds(lBounds.getSize(), 1.0);
  transformPoint(v, lBounds, uBounds, newlBounds, newuBounds);
  return interpolant->eval(v);
}

double SparseGridResponseSurfaceBspline::evalGradientNonUniform(sgpp::base::DataVector v,
                                                                sgpp::base::DataVector& gradient,
                                                                sgpp::base::DataVector lBounds,
                                                                sgpp::base::DataVector uBounds) {
  sgpp::base::DataVector newlBounds(lBounds.getSize(), 0.0);
  sgpp::base::DataVector newuBounds(lBounds.getSize(), 1.0);
  transformPoint(v, lBounds, uBounds, newlBounds, newuBounds);
  return interpolantGradient->eval(v, gradient);
}

double SparseGridResponseSurfaceBspline::getIntegral() {
  sgpp::base::OperationQuadrature* opQuad = sgpp::op_factory::createOperationQuadrature(*grid);
  return opQuad->doQuadrature(coefficients);
}

double SparseGridResponseSurfaceBspline::getMean(sgpp::base::DistributionsVector pdfs,
                                                 size_t quadOrder) {
  sgpp::base::OperationWeightedQuadrature* opWQuad =
      sgpp::op_factory::createOperationWeightedQuadrature(*grid);
  return opWQuad->doWeightedQuadrature(coefficients, pdfs, quadOrder);
}

double SparseGridResponseSurfaceBspline::getVariance(sgpp::base::DistributionsVector pdfs,
                                                     size_t quadOrder) {
  double mean = getMean(pdfs, quadOrder);
  sgpp::base::OperationWeightedSecondMoment* opWSM =
      sgpp::op_factory::createOperationWeightedSecondMoment(*grid);
  double meanSquare = opWSM->doWeightedQuadrature(coefficients, pdfs, quadOrder);
  double variance = meanSquare - mean * mean;
  return variance;
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

  double fX0 = functionValues[x0Index];
  double ftX0 = interpolant->eval(x0);

  std::cout << "x0 = " << x0.toString() << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  sgpp::optimization::optimizer::GradientDescent gradientDescent(*interpolant,
                                                                 *interpolantGradient);
  gradientDescent.setStartingPoint(x0);
  gradientDescent.optimize();
  const sgpp::base::DataVector& xOpt = gradientDescent.getOptimalPoint();
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
    functionValues[i] = objectiveFunc->eval(p);
  }
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "Solving failed!" << std::endl;
  }
}

void SparseGridResponseSurfaceBspline::transformPoint(sgpp::base::DataVector& v,
                                                      sgpp::base::DataVector lBounds,
                                                      sgpp::base::DataVector uBounds,
                                                      sgpp::base::DataVector newlBounds,
                                                      sgpp::base::DataVector newuBounds) {
  v.sub(lBounds);
  uBounds.sub(lBounds);
  v.componentwise_div(uBounds);
  newuBounds.sub(newlBounds);
  v.componentwise_mult(newuBounds);
  v.add(newlBounds);
}

}  // namespace optimization
}  // namespace sgpp
