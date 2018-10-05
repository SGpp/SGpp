// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>

class Rand_double {
 public:
  Rand_double(double low, double high)
      : r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine())) {}

  double operator()() { return r(); }

 private:
  std::function<double()> r;
};

namespace sgpp {
namespace optimization {
void ASMatrixNakBspline::buildRegularInterpolant(size_t level) {
  auto regularGrid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
  sgpp::base::GridStorage& gridStorage = regularGrid->getStorage();
  regularGrid->getGenerator().regular(level);
  evaluationPoints.resizeZero(gridStorage.getSize(), numDim);
  functionValues.resizeZero(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    evaluationPoints.setRow(i, gridPointVector);
    functionValues[i] = objectiveFunc.eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*regularGrid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  this->coefficients = alpha;
  this->grid = regularGrid;
  this->interpolantFlag = 1;
}

void ASMatrixNakBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixNakBspline::createMatrixMonteCarlo(size_t numPoints) {
  if (interpolantFlag == 0) {
    std::cout << "ASMatrixNakBspline: cannot create Matrix without interpolant!\n";
    return;
  }
  sgpp::optimization::InterpolantScalarFunctionGradient interpolantGradient(*(this->grid),
                                                                            this->coefficients);

  //  RandomNumberGenerator::getInstance().setSeed();
  this->C.resize(this->numDim, this->numDim);
  this->C.setZero();
  Rand_double rd{0, 1};
  for (size_t i = 0; i < numPoints; ++i) {
    sgpp::base::DataVector randomVector(this->numDim, 1);
    // todo (rehmemk) somehow the randomnumbergenerator results are much worse than the Rand_double
    // result
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    //    for (size_t d = 0; d < numDim; ++d) {
    //      randomVector[d] = rd();
    //    }
    //    std::cout << randomVector.toString() << std::endl;
    sgpp::base::DataVector gradient(this->numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

double ASMatrixNakBspline::matrixEntryGauss(size_t i, size_t j) {
  double entry = 0.0;
  for (size_t k = 0; k < coefficients.getSize(); k++) {
    for (size_t l = 0; l < coefficients.getSize(); l++) {
      entry += coefficients[k] * coefficients[l] * integralDxbiDxbj(i, j, k, l);
    }
  }
  return entry;
}

double ASMatrixNakBspline::integralDxbiDxbj(size_t i, size_t j, size_t k, size_t l) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double integral = 1.0;
  for (size_t d = 0; d < numDim; d++) {
    sgpp::base::GridPoint& gpk = gridStorage.getPoint(k);
    sgpp::base::GridPoint& gpl = gridStorage.getPoint(l);

    size_t indexk = gpk.getIndex(d);
    size_t indexl = gpl.getIndex(d);
    size_t levelk = gpk.getLevel(d);
    size_t levell = gpl.getLevel(d);

    if ((d == i && i != j)) {
      // int d/dxi b_{k_i} (x_i) b_{l_i} (x_i) dxi
    } else if (d == j && i != j) {
      // int b_{k_j} (x_j) d/dxj b_{l_j} (x_j) dxj
    } else if (d == i && i == j) {
      // int d/dxi b_{k_i} (xi) d/dxi b_{l_i} (xi) dxi
    } else {
      // int b_{k_d} (x_d) b_{l_d} (x_d) dxd
    }
  }
  return 0;
}

double ASMatrixNakBspline::univariateIntegral(size_t level1, size_t index1, size_t index2,
                                              size_t level2, bool dx1, bool dx2) {
  return 0;
}

void ASMatrixNakBspline::createMatrixGauss() {
  C.resize(numDim, numDim);
  for (unsigned int i = 0; i <= C.cols() - 1; i++) {
    for (unsigned int j = i + 1; j < C.rows(); j++) {
      C(i, j) = 1.0;
      C(j, i) = 2.0;
    }
  }
  for (unsigned int i = 0; i < C.rows(); i++) {
    C(i, i) = 3.0;
  }
  std::cout << "C:\n" << C << std::endl;
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
