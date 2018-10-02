// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

/**
 * Class that creates and holds the matrix C for active subspace detection for objective function f
 * C = \int \nabla f \nabla f^T dx
 *
 * Currently only a B-Spline based version exists.
 * It would be desireable to create a Monte Carlo one
 * Also does not yet take probability density functions into account
 * <=> assumes only equally distributed random variables are used
 */
class ASMatrix {
 public:
  /**
   * Constructor
   *
   * @param objectiveFunc the objective function
   */
  ASMatrix(WrapperScalarFunction objectiveFunc) : objectiveFunc(objectiveFunc) {}

  /**
   * Destructor
   */
  virtual ~ASMatrix() {}

  /*
   * Calculates the entries for the matrix C from numPoints samples / a grid with numPoints points
   *
   * @param numPoints number of evaluation points
   */
  virtual void createMatrix(size_t numPoints) = 0;

  /**
   * Uses the Eigen library to calculate an eigenvalue decomposition of C,
   * C = W \Lambda W
   * with \Lambda = diag(\lambda_1, \dots , \lambda_m} the eigenvalues and W the matrix of (column
   * wise) eigen vectors
   */
  void evDecomposition();

  Eigen::VectorXd getEigenvalues() { return this->eigenvalues; };

  /**
   * The Matrix W_1 containing the n first columns of W spans the active subset
   *
   * @param n active subspace indicator (active variables: x_0,\dots , x_{n-1}
   * @return matrix W1
   */
  Eigen::MatrixXd getTransformationMatrix(size_t n) {
    return this->W.block(0, 0, this->W.cols(), n);
  };

 protected:
  /**
   * converts a SG++ DataVector to Eigen vector
   *
   * @param v SG++ DataVector
   * @return Eigen library vector containing the elements of v
   */
  Eigen::VectorXd DataVectorToEigen(sgpp::base::DataVector v) {
    return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
  }

  /**
   * converts an Eigen vector to a SG++ DataVector
   *
   * @param v Eigen library vector
   * @return  SG++ DataVector containing the elements of v
   */
  sgpp::base::DataVector EigenToDataVector(Eigen::VectorXd e) {
    sgpp::base::DataVector v;
    v.resize(e.size());
    Eigen::VectorXd::Map(&v[0], e.size()) = e;
    return v;
  }

  // objective function
  WrapperScalarFunction objectiveFunc;
  // active subspace matrix C = \int \nabla f \nabla f^T dx
  Eigen::MatrixXd C;
  // eigenvectors of C (one per column)
  Eigen::MatrixXd W;
  // eigenvalues of C
  Eigen::VectorXd eigenvalues;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
