// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DENSITYSYSTEMMATRIX_HPP
#define DENSITYSYSTEMMATRIX_HPP

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
class DensitySystemMatrix : public sgpp::base::OperationMatrix {
 private:
  /// the lambda, the regularisation parameter
  double lambda;
  /// Operation A for calculating the data matrix
  /// (L2 Dot-Product of basis functions)
  sgpp::base::OperationMatrix* A;
  /// OperationB for calculating the data matrix
  sgpp::base::OperationMultipleEval* B;
  /// OperationMatrix, the regularisation method
  sgpp::base::OperationMatrix* C;
  /// Training data
  sgpp::base::DataMatrix* data;

 public:
  /**
   * Std-Constructor
   *
   * @param grid  reference to the sparse grid
   * @param trainData reference to DataVector that contains the training data
   * @param C the regression functional
   * @param lambdaRegression the regression parameter
   */
  DensitySystemMatrix(sgpp::base::Grid& grid, sgpp::base::DataMatrix& trainData,
                      sgpp::base::OperationMatrix& C, double lambdaRegression);

  /**
   * Generates the left hand side of the classification equation
   *
   * @param alpha parameters for the sparse grid functions
   * @param result reference to the vector which will contain the result
   */
  void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  /**
   * Generates the right hand side of the classification equation
   *
   * @param b reference to the vector which will contain the result of the
   * matrix vector multiplication on the rhs
   */
  void generateb(sgpp::base::DataVector& b);

  /**
   * Std-Destructor
   */
  virtual ~DensitySystemMatrix();
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DENSITYSYSTEMMATRIX_HPP */
