// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class OperationMatrix for the application of classification for
 * the Systemmatrix by using a density function
 */
class DensitySystemMatrix : public base::OperationMatrix {
 private:
  /// Operation A for calculating the data matrix (L2 Dot-Product of basis functions)
  std::unique_ptr<base::OperationMatrix> A;
  /// OperationB for calculating the data matrix
  std::unique_ptr<base::OperationMultipleEval> B;
  /// OperationMatrix, the regularisation method
  std::unique_ptr<base::OperationMatrix> C;
  /// the lambda, the regularisation parameter
  double lambda;
  /// number of training samples
  size_t numSamples;

 public:
  /**
   * Std-Constructor
   *
   * @param A L^2 dot product matrix of some grid
   * @param B MultipleEval matrix of grid and data points
   * @param C the regression functional
   * @param lambda the regression parameter
   * @param numSamples number of data samples
   */
  DensitySystemMatrix(sgpp::base::OperationMatrix* A, sgpp::base::OperationMultipleEval* B,
                      sgpp::base::OperationMatrix* C, double lambda, size_t numSamples);

  /**
   * Std-Constructor
   *
   * @param grid  reference to the sparse grid
   * @param trainData reference to DataVector that contains the training data
   * @param C the regression functional
   * @param lambda the regression parameter
   */
  DensitySystemMatrix(base::Grid& grid, base::DataMatrix& trainData, base::OperationMatrix* C,
                      double lambda);

  /**
   * Generates the left hand side of the classification equation
   *
   * @param alpha parameters for the sparse grid functions
   * @param result reference to the vector which will contain the result
   */
  void mult(base::DataVector& alpha, base::DataVector& result);

  /**
   * Generates the right hand side of the classification equation
   *
   * @param b reference to the vector which will contain the result of the matrix vector
   * multiplication on the rhs
   */
  void generateb(base::DataVector& b);

  /**
   * Computes the unweighted right hand side of the classification equation
   *
   * @param b reference to the vector which will contain the result of the matrix vector
   * multiplication on the rhs
   */
  void computeUnweightedRhs(base::DataVector& b);

  /**
   * Std-Destructor
   */
  virtual ~DensitySystemMatrix();
};

}  // namespace datadriven
}  // namespace sgpp
