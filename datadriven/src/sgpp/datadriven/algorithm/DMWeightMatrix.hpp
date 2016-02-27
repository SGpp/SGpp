// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMWEIGHTMATRIX_HPP
#define DMWEIGHTMATRIX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix with weight
 */
class DMWeightMatrix : public sgpp::base::OperationMatrix {
 private:
  /// the lambda, the regularisation parameter
  double lamb;
  /// sgpp::base::OperationMatrix, the regularisation mehtod
  sgpp::base::OperationMatrix* C;
  /// OperationB for calculating the data matrix
  sgpp::base::OperationMultipleEval* B;
  /// Pointer to the data vector
  sgpp::base::DataMatrix* data;
  /// Pointer to the weight vector
  sgpp::base::DataVector* weight;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to sgpp::base::DataVector that contains the training data
   * @param C the regression functional
   * @param lambda the lambda, the regression parameter
   * @param w the weights to the training data
   */
  DMWeightMatrix(sgpp::base::Grid& SparseGrid, sgpp::base::DataMatrix& trainData,
                 sgpp::base::OperationMatrix& C, double lambda, sgpp::base::DataVector& w);

  /**
   * Std-Destructor
   */
  virtual ~DMWeightMatrix();

  virtual void mult(sgpp::base::DataVector& alpha,
                    sgpp::base::DataVector& result);

  /**
   * Generates the right hand side of the classification equation
   *
   * @param classes the class information of the training data
   * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
   */
  void generateb(sgpp::base::DataVector& classes, sgpp::base::DataVector& b);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DMWEIGHTMATRIX_HPP */

