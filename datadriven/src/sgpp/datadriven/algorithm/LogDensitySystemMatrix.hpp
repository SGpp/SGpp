// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LOGDENSITYSYSTEMMATRIX_HPP
#define LOGDENSITYSYSTEMMATRIX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix by using a
 * density function
 */
class LogDensitySystemMatrix : public base::OperationMatrix {
 private:
  /// sparse grid reference function
  base::Grid& grid;
  base::DataVector& alpha;
  /// Training data
  base::DataMatrix& data;
  /// the lambda, the regularisation parameter
  float_t lambda;
  /// OperationMatrix, the regularisation method
  base::OperationMatrix& C;
  /// Operation A for calculating the data matrix
  /// (L2 Dot-Product of basis functions)
  std::unique_ptr<base::OperationMatrix> A;
  /// OperationB for calculating the data matrix
  std::unique_ptr<base::OperationMultipleEval> B;

 public:
  /**
   * Std-Constructor
   *
   * @param grid  reference to the sparse grid
   * @param trainData reference to DataVector that contains the training data
   * @param C the regression functional
   * @param lambdaRegression the regression parameter
   */
  LogDensitySystemMatrix(base::Grid& grid, base::DataVector& alphaRef, base::DataMatrix& trainData,
                         base::OperationMatrix& C, double lambdaRegression);

  /**
   * Std-Destructor
   */
  virtual ~LogDensitySystemMatrix();

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
   * @param b reference to the vector which will contain the result of the
   * matrix vector multiplication on the rhs
   */
  void generateb(base::DataVector& b);
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* DENSITYSYSTEMMATRIX_HPP */
