// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/activeSubspaces/EigenFunctionalities.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

/**
 * Class that creates and holds the matrix C for active subspace detection for objective function f
 * C = \int \nabla f \nabla f^T dx
 *
 * Does not yet take probability density functions into account
 * <=> assumes only equally distributed random variables are used
 */
class ASMatrix {
 public:
  /**
   * Constructor
   *
   */
  ASMatrix(sgpp::base::DataMatrix evaluationPoints = sgpp::base::DataMatrix(),
           sgpp::base::DataVector functionValues = sgpp::base::DataVector())
      : evaluationPoints(evaluationPoints), functionValues(functionValues) {}

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

  Eigen::VectorXd getEigenvalues() { return eigenvalues; }
  sgpp::base::DataVector getEigenvaluesDataVector() { return EigenToDataVector(eigenvalues); }

  Eigen::MatrixXd getEigenvectors() { return W; }
  sgpp::base::DataMatrix getEigenvectorsDataMatrix() { return EigenToDataMatrix(W); }

  sgpp::base::DataMatrix getEvaluationPoints() { return evaluationPoints; }

  sgpp::base::DataVector getFunctionValues() { return functionValues; }

  /**
   * The Matrix W_1, containing the n last columns of W, spans the active subset
   * (Eigen sorts the eigenvectors in v increasing)
   *
   * @param n active subspace indicator (active variables: x_n,\dots , x_D
   * @return matrix W1 as Eigen matrix
   */
  Eigen::MatrixXd getTransformationMatrix(size_t n) {
    return W.block(0, W.cols() - n, W.rows(), n);
  }

  /**
   * The Matrix W_1, containing the n last columns of W, spans the active subset
   * (Eigen sorts the eigenvectors in v increasing)
   *
   * @param n active subspace indicator (active variables: x_n,\dots , x_D
   * @return matrix W1 as DataMatrix
   */
  sgpp::base::DataMatrix getTransformationMatrixDataMatrix(size_t n) {
    Eigen::MatrixXd W1 = getTransformationMatrix(n);
    return EigenToDataMatrix(W1);
  }

  Eigen::MatrixXd getMatrix() { return C; }

  sgpp::base::DataMatrix getMatrixDataMatrix() { return EigenToDataMatrix(C); }

  void setMatrix(Eigen::MatrixXd newC) { C = newC; }

 protected:
  // active subspace matrix C = \int \nabla f \nabla f^T dx \in R^{numDim x numDim}
  Eigen::MatrixXd C;
  // eigenvectors of C (one per column)
  Eigen::MatrixXd W;
  // eigenvalues of C
  Eigen::VectorXd eigenvalues;
  // points the objective function is evaluated at (row-wise)
  sgpp::base::DataMatrix evaluationPoints;
  // evaluations of the objective function
  sgpp::base::DataVector functionValues;
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
