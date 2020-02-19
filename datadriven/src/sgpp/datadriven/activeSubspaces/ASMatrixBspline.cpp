// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/datadriven/activeSubspaces/ASMatrixBspline.hpp>

namespace sgpp {
namespace datadriven {

void ASMatrixBspline::initialize(sgpp::base::GridType gridType) {
  this->gridType = gridType;
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::ModNakBspline) {
    grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

void ASMatrixBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixBspline::createMatrixMonteCarlo(size_t numMCPoints) {
  sgpp::base::InterpolantScalarFunctionGradient interpolantGradient(*grid, coefficients);

  sgpp::base::RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  for (size_t i = 0; i < numMCPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    sgpp::base::RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    sgpp::base::DataVector gradient(numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numMCPoints);
}

void ASMatrixBspline::createMatrixGauss() {
  C.resize(numDim, numDim);
  {
    size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0)) * 2;
    sgpp::datadriven::NakBsplineScalarProducts scalarProducts(gridType, gridType, degree, degree,
                                                              quadOrder);
    for (unsigned int i = 0; i <= C.cols(); i++) {
      for (unsigned int j = i; j < C.rows(); j++) {
        double entry = matrixEntryGauss(i, j, scalarProducts);
        C(i, j) = entry;
        C(j, i) = entry;
      }
    }
  }
}

void ASMatrixBspline::toFile(std::string path) {
  std::string gridPath = path + "/ASMGrid.grid";
  std::filebuf fb;
  fb.open(gridPath, std::ios::out);
  std::ostream os(&fb);
  os << grid->serialize();
  fb.close();
  std::string coeffPath = path + "/ASMCoefficients.vec";
  coefficients.toFile(coeffPath);
}

sgpp::datadriven::ASResponseSurfaceNakBspline ASMatrixBspline::getResponseSurfaceInstance(
    size_t asDimension, sgpp::base::GridType gridType, size_t degree) {
  Eigen::MatrixXd W1 = this->getTransformationMatrix(asDimension);

  sgpp::datadriven::ASResponseSurfaceNakBspline responseSurf(W1, gridType, degree);
  return responseSurf;
}

// ----------------- auxiliary routines -----------

void ASMatrixBspline::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  this->calculateCoefficients();
}

double ASMatrixBspline::matrixEntryGauss(
    size_t i, size_t j, sgpp::datadriven::NakBsplineScalarProducts& scalarProducts) {
  double entry = 0.0;
  // ToDo (rehmemk) Does the parallelization work fine? If yes remove the scalarProducts and use
  // only the innerScalarProducts
  size_t numBasisFunctions = coefficients.getSize();
#pragma omp parallel
  {
    //    std::cout << "ASM using" << omp_get_num_threads() << " threads\n";
    size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0)) * 2;
    sgpp::datadriven::NakBsplineScalarProducts innerScalarProducts(gridType, gridType, degree,
                                                                   degree, quadOrder);
#pragma omp for reduction(+ : entry) collapse(2)
    for (size_t k = 0; k < numBasisFunctions; k++) {
      for (size_t l = 0; l < numBasisFunctions; l++) {
        entry += coefficients[k] * coefficients[l] *
                 scalarProductDxbiDxbj(i, j, k, l, innerScalarProducts);
      }
    }
  }
  return entry;
}

double ASMatrixBspline::scalarProductDxbiDxbj(
    size_t i, size_t j, size_t k, size_t l,
    sgpp::datadriven::NakBsplineScalarProducts& scalarProducts) {
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
      integral1D = scalarProducts.basisScalarProduct(levelk, indexk, true, levell, indexl, false);
    } else if (d == j && i != j) {
      // int b_{k_j} (x_j) d/dxj b_{l_j} (x_j) dxj
      integral1D = scalarProducts.basisScalarProduct(levelk, indexk, false, levell, indexl, true);
    } else if (d == i && i == j) {
      // int d/dxi b_{k_i} (xi) d/dxi b_{l_i} (xi) dxi
      integral1D = scalarProducts.basisScalarProduct(levelk, indexk, true, levell, indexl, true);
    } else {
      // int b_{k_d} (x_d) b_{l_d} (x_d) dxd
      integral1D = scalarProducts.basisScalarProduct(levelk, indexk, false, levell, indexl, false);
    }
    if (integral1D == 0) return 0.0;
    integral *= integral1D;
  }
  return integral;
}

double ASMatrixBspline::evalInterpolant(sgpp::base::DataVector v) {
  sgpp::base::InterpolantScalarFunction interpolant(*grid, coefficients);
  return interpolant.eval(v);
}

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
