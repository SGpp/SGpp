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
 * Class that implements the virtual class OperationMatrix for the application of density estimation
 * for the Systemmatrix by using a density difference function for two input datasets
 */
class DensityDifferenceSystemMatrix : public base::OperationMatrix {
 private:
  /// Operation A for calculating the data matrix
  /// (L2 Dot-Product of basis functions)
  std::unique_ptr<base::OperationMatrix> A;
  /// OperationB for calculating the data matrix for first dataset
  std::unique_ptr<base::OperationMultipleEval> B_p;
  /// OperationB for calculating the data matrix for second dataset
  std::unique_ptr<base::OperationMultipleEval> B_q;
  /// OperationMatrix, the regularisation method
  std::unique_ptr<base::OperationMatrix> C;
  /// the lambda, the regularisation parameter
  double lambda;
  /// number of training samples
  size_t numSamples_P, numSamples_Q;

 public:
  /**
   * Std-Constructor
   *
   * @param A L^2 dot product matrix of some grid
   * @param B_p MultipleEval matrix of grid and data points for first dataset
   * @param B_q MultipleEval matrix of grid and data points for second dataset
   * @param C the regression functional
   * @param lambda the regression parameter
   * @param numSamples_P number of data samples for first dataset
   * @param numSamples_Q number of data samples for second dataset
   */
  DensityDifferenceSystemMatrix(sgpp::base::OperationMatrix* A,
                                sgpp::base::OperationMultipleEval* B_p,
                                sgpp::base::OperationMultipleEval* B_q,
                                sgpp::base::OperationMatrix* C, double lambda, size_t numSamples_P,
                                size_t numSamples_Q);

  /**
   * Std-Constructor
   *
   * @param grid  reference to the sparse grid
   * @param trainData_P reference to DataVector that contains the training data from first dataset
   * @param trainData_Q reference to DataVector that contains the training data from second dataset
   * @param C the regression functional
   * @param lambda the regression parameter
   */
  DensityDifferenceSystemMatrix(base::Grid& grid, base::DataMatrix& trainData_P,
                                base::DataMatrix& trainData_Q, base::OperationMatrix* C,
                                double lambda);

  /**
   * Generates the left hand side of the density estimation equation
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
   * Computes the unweighted right hand sides of the density estimation equation
   *
   * @param bp reference to the vector which will contain the result of the matrix vector
   * multiplication on the rhs for first dataset
   * @param bq reference to the vector which will contain the result of the matrix vector
   * multiplication on the rhs for second dataset
   */
  void computeUnweightedRhs(base::DataVector& bp, base::DataVector& bq);

  /**
   * Std-Destructor
   */
  virtual ~DensityDifferenceSystemMatrix();
};

}  // namespace datadriven
}  // namespace sgpp
