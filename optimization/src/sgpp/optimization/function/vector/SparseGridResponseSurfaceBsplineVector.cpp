// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>

#include <sgpp/optimization/function/vector/SparseGridResponseSurfaceBsplineVector.hpp>

namespace sgpp {
namespace optimization {

void SparseGridResponseSurfaceBsplineVector::regular(size_t level) {
  grid->getGenerator().regular(level);
  calculateInterpolationCoefficients();
  interpolants = std::make_unique<sgpp::base::InterpolantVectorFunction>(*grid, coefficients);
  interpolantGradients =
      std::make_unique<sgpp::base::InterpolantVectorFunctionGradient>(*grid, coefficients);
}

void SparseGridResponseSurfaceBsplineVector::regularByPoints(size_t numPoints, bool verbose) {
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
  interpolants = std::make_unique<sgpp::base::InterpolantVectorFunction>(*grid, coefficients);
  interpolantGradients =
      std::make_unique<sgpp::base::InterpolantVectorFunctionGradient>(*grid, coefficients);
}

void SparseGridResponseSurfaceBsplineVector::surplusAdaptive(size_t maxNumGridPoints,
                                                             size_t initialLevel,
                                                             size_t refinementsNum, bool verbose) {
  regular(initialLevel);
  while (grid->getSize() < maxNumGridPoints) {
    refineSurplusAdaptive(refinementsNum);
    if (verbose)
      std::cout << "Refining. Calculated a grid with " << grid->getSize() << " points.\n";
  }
  interpolants = std::make_unique<sgpp::base::InterpolantVectorFunction>(*grid, coefficients);
  interpolantGradients =
      std::make_unique<sgpp::base::InterpolantVectorFunctionGradient>(*grid, coefficients);
}

// void SparseGridResponseSurfaceBsplineVector::ritterNovak(size_t maxNumGridPoints, double gamma,
//                                                          bool verbose) {
//   // this uses the default values for Ritter Novaks initial Level, max level, powerMethod
//   sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(*objectiveFunc, *grid,
//                                                                 maxNumGridPoints, gamma);
//   if (!gridGen.generate()) {
//     std::cout << "Grid generation failed, exiting.\n";
//   }
//   calculateInterpolationCoefficients();
//   interpolant =
//       std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
//   interpolantGradient =
//   std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
//       *grid, coefficients);
//   if (verbose) std::cout << "Refining. Calculated a grid with " << grid->getSize() << "
//   points.\n";
// }

sgpp::base::DataVector SparseGridResponseSurfaceBsplineVector::eval(sgpp::base::DataVector v) {
  transformPoint(v, lb, ub, unitLBounds, unitUBounds);
  sgpp::base::DataVector evaluations(numRes);
  interpolants->eval(v, evaluations);
  return evaluations;
}
sgpp::base::DataVector SparseGridResponseSurfaceBsplineVector::evalJacobian(
    sgpp::base::DataVector v, sgpp::base::DataMatrix& jacobian) {
  jacobian.resizeZero(numRes, numDim);
  transformPoint(v, lb, ub, unitLBounds, unitUBounds);
  sgpp::base::DataVector evaluations(numRes);
  interpolantGradients->eval(v, evaluations, jacobian);
  // scale gradient components according to the inner derivative of the chain rule when transforming
  // the interpolation point from the original coordinates to unit cube

  for (size_t i = 0; i < numRes; i++) {
    for (size_t j = 0; j < numDim; j++) {
      double temp = jacobian.get(i, j) / (ub[j] - lb[j]);
      jacobian.set(i, j, temp);
    }
  }
  return evaluations;
}

sgpp::base::DataVector SparseGridResponseSurfaceBsplineVector::getIntegrals() {
  sgpp::base::OperationQuadrature* opQuad = sgpp::op_factory::createOperationQuadrature(*grid);
  sgpp::base::DataVector integrals(numRes);
  double vol = domainVolume();
  for (size_t t = 0; t < numRes; t++) {
    sgpp::base::DataVector coefficients_t(grid->getSize());
    coefficients.getColumn(t, coefficients_t);
    integrals[t] = opQuad->doQuadrature(coefficients_t);
    integrals[t] *= vol;
  }
  return integrals;
}

// double SparseGridResponseSurfaceBsplineVector::getMeans(sgpp::base::DistributionsVector pdfs,
//                                                        size_t quadOrder) {
//   sgpp::base::OperationWeightedQuadrature* opWQuad =
//       sgpp::op_factory::createOperationWeightedQuadrature(*grid);
//   mean = opWQuad->doWeightedQuadrature(coefficients, pdfs, quadOrder);
//   computedMeanFlag = true;
//   return mean;
// }

// sgpp::base::DataVector SparseGridResponseSurfaceBsplineVector::getVariances(
//     sgpp::base::DistributionsVector pdfs, size_t quadOrder) {
//   if (!computedMeanFlag) {
//     double dummy = getMean(pdfs, quadOrder);
//   }
//   sgpp::base::OperationWeightedSecondMoment* opWSM =
//       sgpp::op_factory::createOperationWeightedSecondMoment(*grid);
//   double meanSquare = opWSM->doWeightedQuadrature(coefficients, pdfs, quadOrder);
//   std::cout << std::setprecision(16) << "Variance Calculation. Mean: " << mean
//             << " meanSquare: " << meanSquare << "\n";
//   double variance = meanSquare - mean * mean;
//   sgpp::base::DataVector returnVec(3);
//   returnVec[0] = variance;
//   returnVec[1] = meanSquare;
//   returnVec[2] = mean;
//   return returnVec;
// }

// sgpp::base::DataVector SparseGridResponseSurfaceBsplineVector::optimize() {
//   /**
//    * The gradient method needs a starting point.
//    *  We use the grid point with the smallest (most promising) function value and save it in x0.
//    */
//   sgpp::base::DataVector x0(numDim);
//   sgpp::base::GridStorage& gridStorage = grid->getStorage();

//   // index of grid point with minimal function value
//   size_t x0Index =
//       std::distance(functionValues.getPointer(),
//                     std::min_element(functionValues.getPointer(),
//                                      functionValues.getPointer() + functionValues.getSize()));

//   x0 = gridStorage.getCoordinates(gridStorage[x0Index]);

//   sgpp::optimization::optimizer::GradientDescent gradientDescent(*interpolant,
//                                                                  *interpolantGradient);
//   gradientDescent.setStartingPoint(x0);
//   gradientDescent.optimize();
//   const sgpp::base::DataVector& unitXOpt = gradientDescent.getOptimalPoint();
//   sgpp::base::DataVector xOpt(unitXOpt);
//   transformPoint(xOpt, unitLBounds, unitUBounds, lb, ub);
//   return xOpt;
// }

std::string SparseGridResponseSurfaceBsplineVector::serializeGrid() {
  std::string gridStr;
  grid->serialize(gridStr);
  return gridStr;
}

// ----------------- auxiliary routines -----------

void SparseGridResponseSurfaceBsplineVector::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::VectorSurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  calculateInterpolationCoefficients();
}

void SparseGridResponseSurfaceBsplineVector::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  functionValues.resizeZero(gridStorage.getSize(), numRes);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::DataVector p = gridStorage.getPointCoordinates(i);
    transformPoint(p, unitLBounds, unitUBounds, lb, ub);
    sgpp::base::DataVector evaluations(numRes);
    objectiveFunc->eval(p, evaluations);
    functionValues.setRow(i, evaluations);
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
