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
#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/Node.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix by using a
 * density function
 */
class PiecewiseConstantSmoothedRegressionSystemMatrix: public
  sgpp::base::OperationMatrix {
 private:
  sgpp::datadriven::PiecewiseConstantRegression::Node& piecewiseRegressor;
  sgpp::base::Grid& grid;
  /// the lambda, the regularisation parameter
  double lambda;
  /// Operation A for calculating the data matrix
  /// (L2 Dot-Product of basis functions)
  sgpp::base::OperationMatrix* A;
  /// OperationB for calculating the data matrix
  //        sgpp::base::OperationMultipleEval* B;
  /// OperationMatrix, the regularisation method
  sgpp::base::OperationMatrix* C;

 public:
  /**
   * Std-Constructor
   *
   * @param piecewiseRegressor approximation with piecewise-constant octtree
   * @param grid  reference to the sparse grid
   * @param C the regression functional
   * @param lambdaRegression the regression parameter
   */
  PiecewiseConstantSmoothedRegressionSystemMatrix(
    sgpp::datadriven::PiecewiseConstantRegression::Node& piecewiseRegressor,
    sgpp::base::Grid& grid,
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
  virtual ~PiecewiseConstantSmoothedRegressionSystemMatrix();
};

}  // namespace datadriven
}  // namespace sgpp

