// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/activeSubspaces/EigenFunctionalities.hpp>
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
  ASMatrix(WrapperScalarFunction objectiveFunc)
      : objectiveFunc(objectiveFunc), numDim(objectiveFunc.getNumberOfParameters()) {}

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
  void evDecompositionForSymmetricMatrices();

  Eigen::VectorXd getEigenvalues() { return eigenvalues; };

  Eigen::MatrixXd getEigenvectors() { return W; };

  sgpp::base::DataMatrix getEvaluationPoints() { return evaluationPoints; }

  sgpp::base::DataVector getFunctionValues() { return functionValues; }

  /**
   * The Matrix W_1 containing the n last columns of W spans the active subset
   *
   * @param n active subspace indicator (active variables: x_n,\dots , x_D
   * @return matrix W1
   */
  Eigen::MatrixXd getTransformationMatrix(size_t n) { return W.block(0, n, W.cols(), numDim - n); };

  Eigen::MatrixXd getMatrix() { return C; }

  void setMatrix(Eigen::MatrixXd newC);

 protected:
  // objective function
  WrapperScalarFunction objectiveFunc;
  // dimensionality
  size_t numDim;
  // active subspace matrix C = \int \nabla f \nabla f^T dx \in R^{numDim x numDim}
  Eigen::MatrixXd C;
  // eigenvectors of C (one per column)
  Eigen::MatrixXd W;
  // eigenvalues of C
  Eigen::VectorXd eigenvalues;
  // points the objective function is evaluated at (each row is one point)
  sgpp::base::DataMatrix evaluationPoints;
  // evaluations of the objective function
  sgpp::base::DataVector functionValues;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
