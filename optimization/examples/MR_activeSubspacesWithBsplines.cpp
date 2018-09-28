// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>

#include <fstream>
#include <functional>
#include <iostream>
#include <random>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

class objectiveFunction {
 public:
  objectiveFunction() {
    alpha = 1000;
    beta = 1;
  }

  double eval(sgpp::base::DataVector x) {
    double res = alpha * (x[0] * x[0] + x[2] * x[2]) + beta * (x[1] * x[1] + x[3] * x[3]);
    //    double res = exp(0.7 * x[0] + 0.3 * x[1]);
    return res;
  }

  sgpp::base::DataVector grad(sgpp::base::DataVector x) {
    sgpp::base::DataVector res(x.getSize());
    //    res[0] = 2 * alpha * x[0];
    //    res[1] = 2 * beta * x[1];
    //    res[2] = 2 * alpha * x[2];
    //    res[3] = 2 * beta * x[3];
    res[0] = 0.7 * exp(0.7 * x[0] + 0.3 * x[1]);
    res[1] = 0.3 * exp(0.7 * x[0] + 0.3 * x[1]);
    return res;
  }

 private:
  double alpha;
  double beta;
};

Eigen::VectorXd DataVectorToEigen(sgpp::base::DataVector v) {
  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
}

sgpp::base::DataVector EigenToDataVector(Eigen::VectorXd v1) {
  sgpp::base::DataVector v2;
  v2.resize(v1.size());
  Eigen::VectorXd::Map(&v2[0], v1.size()) = v1;
  return v2;
}
class reducedInterpolant {
 public:
  reducedInterpolant(std::shared_ptr<sgpp::base::Grid> grid, Eigen::MatrixXd W1,
                     sgpp::base::DataVector coefficients)
      : grid(grid), I(*grid, coefficients), W1T(W1.transpose()), coefficients(coefficients) {}

  double eval(sgpp::base::DataVector v) {
    Eigen::VectorXd v_Eigen = DataVectorToEigen(v);
    Eigen::VectorXd trans_v_Eigen = W1T * v_Eigen;
    sgpp::base::DataVector trans_v_DataVector = EigenToDataVector(trans_v_Eigen);
    return I.eval(trans_v_DataVector);
  }

  size_t getDim() { return grid->getDimension(); }

 private:
  std::shared_ptr<sgpp::base::Grid> grid;
  sgpp::optimization::InterpolantScalarFunction I;
  Eigen::MatrixXd W1T;
  sgpp::base::DataVector coefficients;
};

// has no seed. Same random numbers every time
class Rand_double {
 public:
  Rand_double(double low, double high)
      : r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine())) {}

  double operator()() { return r(); }

 private:
  std::function<double()> r;
};

/**
 * return the indices according to a creasing ordering of the inputs
 *
 * @param v vector containing data
 * @return vector of indices according to the decreasing sorting of the input
 */
std::vector<size_t> sortingIndices(Eigen::VectorXd& v) {
  // initialize original index locations
  std::vector<size_t> idx(v.rows());
  std::iota(idx.begin(), idx.end(), 0);
  // sort indices based on comparing values in v
  sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v(i1) > v(i2); });
  return idx;
}

void sortEigenValuesAndVectors(Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
  std::vector<size_t> idx = sortingIndices(eigenvalues);
  Eigen::VectorXd tempEigenvalues(eigenvalues.rows());
  Eigen::MatrixXd tempEigenvectors(eigenvectors.cols(), eigenvectors.rows());
  size_t k = 0;
  for (size_t& i : idx) {
    tempEigenvalues(k) = eigenvalues(i);
    tempEigenvectors.block(0, k, eigenvectors.rows(), 1) =
        eigenvectors.block(0, i, eigenvectors.rows(), 1);
    k++;
  }
  eigenvalues = tempEigenvalues;
  eigenvectors = tempEigenvectors;
}

int main() {
  size_t numDim = 4;
  size_t degree = 3;
  size_t level = 5;
  objectiveFunction objectiveFunc;
  // specify active subspace
  // active variables: x_0,...,x_{n-1} inactive variables x_{n},...,x_m
  size_t n = 1;
  // set k to one more than the desired reduced dimensionality
  size_t k = n + 1;
  // choose an oversampling factor alpha between 2 and 10 (the larger, the more evaluations, the
  // better the approximation of the eigenvalues)
  size_t alpha = 3;
  // the number of samples ensuring that the first k eigenvalues are accurate enough
  size_t M = alpha * k * static_cast<size_t>(ceil(log(static_cast<double>(numDim))));

  // build sparse grid B-spline interpolant
  std::unique_ptr<sgpp::base::Grid> detectionGrid(
      sgpp::base::Grid::createNakBsplineGrid(numDim, degree));
  sgpp::base::GridStorage& detectionGridStorage = detectionGrid->getStorage();
  detectionGrid->getGenerator().regular(level);
  sgpp::base::DataVector functionValues(detectionGridStorage.getSize());
  for (size_t i = 0; i < detectionGridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = detectionGridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(detectionGridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    functionValues[i] = objectiveFunc.eval(gridPointVector);
  }
  sgpp::base::DataVector coeffs(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*detectionGrid);

  sgpp::optimization::sle_solver::Armadillo sleSolver;

  // solve linear system
  if (!sleSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  sgpp::optimization::InterpolantScalarFunction detectionInterpolant(*detectionGrid, coeffs);
  sgpp::optimization::InterpolantScalarFunctionGradient detectionInterpolantGradient(*detectionGrid,
                                                                                     coeffs);

  // This is a Monte carlo approxiation of the integral. Do this analytically exact
  // Precalculate relevant values for speed up
  // Parallelize
  Rand_double rd{0, 1};
  sgpp::base::DataVector randomVector(numDim, 1);
  Eigen::MatrixXd C(numDim, numDim);
  C.setZero();
  for (size_t i = 0; i < M; ++i) {
    for (size_t d = 0; d < numDim; ++d) {
      randomVector[d] = rd();
    }
    //    sgpp::base::DataVector gradient = objectiveFunc.grad(randomVector);
    sgpp::base::DataVector gradient(numDim);
    detectionInterpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    C += e * e.transpose();
  }
  C /= static_cast<double>(M);

  // calculate EW and EV
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
  if (eigensolver.info() != Eigen::Success) abort();
  Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
  Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
  std::cout << "EW:\n" << eigenvalues << std::endl;
  std::cout << "EV(columns):\n" << eigenvectors << std::endl;

  // is this necessary or are eigenvalues always in ascending order?
  sortEigenValuesAndVectors(eigenvalues, eigenvectors);
  std::cout << "EW sorted:\n" << eigenvalues << std::endl;
  std::cout << "EV sorted:\n" << eigenvectors << std::endl;

  //  // write EW and EV to file
  //  std::fstream dataFile;
  //  dataFile.open("/home/rehmemk/git/activeSubspaceSGpp/activeSubSpaces_Python/data/eigenvalues.dat",
  //                std::ofstream::out | std::ofstream::trunc);
  //  for (unsigned int i = 0; i < eigenvalues.rows(); i++) {
  //    dataFile << eigenvalues[i] << "\n";
  //  }
  //  dataFile.close();
  //  dataFile.open("/home/rehmemk/git/activeSubspaceSGpp/activeSubSpaces_Python/data/eigenvectors.dat",
  //                std::ofstream::out | std::ofstream::trunc);
  //  for (unsigned int i = 0; i < eigenvectors.rows(); i++) {
  //    for (unsigned int j = 0; j < eigenvectors.cols() - 1; j++) {
  //      dataFile << eigenvectors(i, j) << ", ";
  //    }
  //    dataFile << eigenvectors(i, eigenvectors.cols() - 1) << "\n";
  //  }
  //  dataFile.close();

  Eigen::MatrixXd W1 = eigenvectors.block(0, 0, eigenvectors.cols(), n);
  Eigen::MatrixXd W2 = eigenvectors.block(0, n, eigenvectors.cols(), numDim - n);

  std::cout << "W1:\n" << W1 << std::endl;
  std::cout << "W2:\n" << W2 << std::endl;

  // Use actual values for debugging
  //  W1(0) = 0.7;
  //  W1(1) = 0.3;

  // create dimension reduced interpolant
  // use grid points of original iterpolant used for active subspace detection
  size_t reducedGridLevel = 5;
  auto reducedGrid = std::make_shared<sgpp::base::NakBsplineGrid>(n, degree);
  sgpp::base::SNakBsplineBase reducedBasis(degree);
  sgpp::base::GridStorage& reducedGridStorage = reducedGrid->getStorage();
  reducedGrid->getGenerator().regular(reducedGridLevel);

  std::cout << "building reduced Matrix" << std::endl;
  Eigen::MatrixXd reducedInterpolationMatrix(detectionGridStorage.getSize(),
                                             reducedGridStorage.getSize());
  for (size_t i = 0; i < detectionGridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gpPoint = detectionGridStorage.getPoint(i);
    // eval in y = W1T * gpPoint
    Eigen::VectorXd gpPoint_trans_Eigen(gpPoint.getDimension());
    for (size_t d = 0; d < gpPoint.getDimension(); d++) {
      gpPoint_trans_Eigen(d) = detectionGridStorage.getUnitCoordinate(gpPoint, d);
    }
    gpPoint_trans_Eigen = W1.transpose() * gpPoint_trans_Eigen;
    for (size_t j = 0; j < reducedGridStorage.getSize(); j++) {
      sgpp::base::GridPoint& gpreducedBasis = reducedGridStorage.getPoint(j);
      double basisEval = 1;
      for (size_t t = 0; t < reducedGridStorage.getDimension(); t++) {
        double basisEval1D = reducedBasis.eval(gpreducedBasis.getLevel(t),
                                               gpreducedBasis.getIndex(t), gpPoint_trans_Eigen(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      reducedInterpolationMatrix(i, j) = basisEval;
    }
  }
  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);
  // Least Squares Fit
  Eigen::VectorXd reducedCoefficients_Eigen =
      reducedInterpolationMatrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(functionValues_Eigen);

  sgpp::base::DataVector reducedCoefficients = EigenToDataVector(reducedCoefficients_Eigen);
  reducedInterpolant rI(reducedGrid, W1, reducedCoefficients);

  sgpp::base::DataVector x(numDim, 0.5);
  std::cout << objectiveFunc.eval(x) << std::endl;
  std::cout << detectionInterpolant.eval(x) << std::endl;
  std::cout << rI.eval(x) << std::endl;

  return 0;
}
