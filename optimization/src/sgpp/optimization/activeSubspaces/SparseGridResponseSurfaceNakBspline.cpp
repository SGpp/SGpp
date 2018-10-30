// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/SparseGridResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace optimization {

void SparseGridResponseSurfaceNakBspline::initialize() {
  numDim = objectiveFunc->getNumberOfParameters();
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_unique<sgpp::base::NakBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_unique<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = std::make_unique<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

void SparseGridResponseSurfaceNakBspline::createRegularResponseSurface(size_t level) {
  grid->getGenerator().regular(level);
  sgpp::base::DataVector alpha = calculateInterpolationCoefficients();
  interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
}

double SparseGridResponseSurfaceNakBspline::eval(sgpp::base::DataVector v) {
  return interpolant->eval(v);
}
double SparseGridResponseSurfaceNakBspline::evalGradient(sgpp::base::DataVector v,
                                                         sgpp::base::DataVector& gradient) {
  return interpolantGradient->eval(v, gradient);
}

void SparseGridResponseSurfaceNakBspline::createSurplusAdaptiveResponseSurface(
    size_t maxNumGridPoints, size_t initialLevel) {
  // number of points to be refined in each step
  size_t refinementsNum = 3;
  grid->getGenerator().regular(initialLevel);
  sgpp::base::DataVector alpha = calculateInterpolationCoefficients();
  interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum, alpha);
    interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
    interpolantGradient =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
  }
}

// ----------------- auxiliary routines -----------

void SparseGridResponseSurfaceNakBspline::refineSurplusAdaptive(size_t refinementsNum,
                                                                sgpp::base::DataVector& alpha) {
  sgpp::base::SurplusRefinementFunctor functor(alpha, refinementsNum);
  grid->getGenerator().refine(functor);
  alpha = calculateInterpolationCoefficients();
}

sgpp::base::DataVector SparseGridResponseSurfaceNakBspline::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    f_values[i] = objectiveFunc->eval(p);
  }
  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  return alpha;
}

}  // namespace optimization
}  // namespace sgpp
