// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>

namespace sgpp {
namespace optimization {
void ASMatrixNakBspline::buildRegularInterpolant(size_t level) {
  grid->getGenerator().regular(level);
  this->calculateInterpolationCoefficients();
}

void ASMatrixNakBspline::buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                                  size_t refinementsNum) {
  grid->getGenerator().regular(initialLevel);
  this->calculateInterpolationCoefficients();
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum);
  }
}

void ASMatrixNakBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixNakBspline::createMatrixMonteCarlo(size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunctionGradient interpolantGradient(*grid, coefficients);

  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  for (size_t i = 0; i < numMCPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    sgpp::base::DataVector gradient(numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numMCPoints);
}

void ASMatrixNakBspline::createMatrixGauss() {
  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0)) * 2;
  sgpp::optimization::NakBsplineScalarProducts scalarProducts(gridType, gridType, degree, degree,
                                                              quadOrder);
  C.resize(numDim, numDim);
  for (unsigned int i = 0; i <= C.cols(); i++) {
    for (unsigned int j = i; j < C.rows(); j++) {
      double entry = matrixEntryGauss(i, j, scalarProducts);
      C(i, j) = entry;
      C(j, i) = entry;
    }
  }
}

double ASMatrixNakBspline::l2InterpolationError(size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunction interpolant(*grid, coefficients);
  double l2Err = 0.0;
  sgpp::base::DataVector randomVector(objectiveFunc->getNumberOfParameters());
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    double evalInterpolant = interpolant.eval(randomVector);
    double evalObjectiveFunc = objectiveFunc->eval(randomVector);
    l2Err += std::pow(evalInterpolant - evalObjectiveFunc, 2.0);
  }
  l2Err = sqrt(l2Err / static_cast<double>(numMCPoints));
  return l2Err;
}

sgpp::base::DataVector ASMatrixNakBspline::l2InterpolationGradientError(
    std::shared_ptr<sgpp::optimization::WrapperScalarFunctionGradient> objectiveFuncGradient,
    size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunctionGradient interpolant(*grid, coefficients);
  size_t numDim = objectiveFuncGradient->getNumberOfParameters();
  sgpp::base::DataVector errors(numDim, 0);
  sgpp::base::DataVector randomVector(numDim);
  sgpp::base::DataVector interpolantEval(numDim);
  sgpp::base::DataVector gradientEval(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    sgpp::optimization::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    interpolant.eval(randomVector, interpolantEval);
    objectiveFuncGradient->eval(randomVector, gradientEval);
    for (size_t d = 0; d < numDim; d++) {
      errors[d] += std::pow(interpolantEval[d] - gradientEval[d], 2.0);
    }
  }
  for (size_t d = 0; d < numDim; d++) {
    errors[d] = sqrt(errors[d] / static_cast<double>(numMCPoints));
  }
  return errors;
}

void ASMatrixNakBspline::toFile(std::string path) {
  std::string gridPath = path + "/ASMGrid.grid";
  std::filebuf fb;
  fb.open(gridPath, std::ios::out);
  std::ostream os(&fb);
  os << grid->serialize();
  fb.close();
  std::string coeffPath = path + "/ASMCoefficients.vec";
  coefficients.toFile(coeffPath);
}

sgpp::optimization::ASResponseSurfaceNakBspline ASMatrixNakBspline::getResponseSurfaceInstance(
    size_t asDimension, sgpp::base::GridType gridType, size_t degree) {
  Eigen::MatrixXd W1 = this->getTransformationMatrix(asDimension);

  sgpp::optimization::ASResponseSurfaceNakBspline responseSurf(W1, gridType, degree);
  return responseSurf;
}

// ----------------- auxiliary routines -----------

void ASMatrixNakBspline::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  this->calculateInterpolationCoefficients();
}

void ASMatrixNakBspline::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  evaluationPoints.resizeZero(gridStorage.getSize(), numDim);
  // ToDo (rehmemk) when refining iteratively use the function values from the last steps!
  functionValues.resizeZero(gridStorage.getSize());
  sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      gridPointVector[d] = gridStorage.getPointCoordinate(i, d);
    }
    evaluationPoints.setRow(i, gridPointVector);
    functionValues[i] = objectiveFunc->eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  coefficients = alpha;
}

double ASMatrixNakBspline::matrixEntryGauss(
    size_t i, size_t j, sgpp::optimization::NakBsplineScalarProducts& scalarProducts) {
  double entry = 0.0;
  for (size_t k = 0; k < coefficients.getSize(); k++) {
    for (size_t l = 0; l < coefficients.getSize(); l++) {
      entry +=
          coefficients[k] * coefficients[l] * scalarProductDxbiDxbj(i, j, k, l, scalarProducts);
    }
  }
  return entry;
}

// Todo(rehmemk) Check already here if (multidimensional) B-spline supports overlap and return 0 if
// so
double ASMatrixNakBspline::scalarProductDxbiDxbj(
    size_t i, size_t j, size_t k, size_t l,
    sgpp::optimization::NakBsplineScalarProducts& scalarProducts) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double integral = 1.0;
  for (size_t d = 0; d < numDim; d++) {
    unsigned int indexk = gridStorage.getPointIndex(k, d);
    unsigned int indexl = gridStorage.getPointIndex(l, d);
    unsigned int levelk = gridStorage.getPointLevel(k, d);
    unsigned int levell = gridStorage.getPointLevel(l, d);
    double integral1D = 0.0;

    if ((d == i && i != j)) {
      // int d/dxi b_{k_i} (x_i) b_{l_i} (x_i) dxi
      integral1D =
          scalarProducts.univariateScalarProduct(levelk, indexk, true, levell, indexl, false);
    } else if (d == j && i != j) {
      // int b_{k_j} (x_j) d/dxj b_{l_j} (x_j) dxj
      integral1D =
          scalarProducts.univariateScalarProduct(levelk, indexk, false, levell, indexl, true);
    } else if (d == i && i == j) {
      // int d/dxi b_{k_i} (xi) d/dxi b_{l_i} (xi) dxi
      integral1D =
          scalarProducts.univariateScalarProduct(levelk, indexk, true, levell, indexl, true);
    } else {
      // int b_{k_d} (x_d) b_{l_d} (x_d) dxd
      integral1D =
          scalarProducts.univariateScalarProduct(levelk, indexk, false, levell, indexl, false);
    }
    if (integral1D == 0) return 0.0;
    integral *= integral1D;
  }
  return integral;
}

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
