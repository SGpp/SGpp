// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXDRE_HPP
#define DMSYSTEMMATRIXDRE_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class base::OperationMatrix for the
 * application of classification for the Systemmatrix
 */
class DMSystemMatrixDRE : public DMSystemMatrixBase {
 private:
  base::Grid& grid;
  /// the rhs dataset
  base::DataMatrix rhs_dataset_;
  /// base::OperationMatrix, the regularisation method
  std::shared_ptr<base::OperationMatrix> C;
  /// OperationB for calculating the data matrix
  std::unique_ptr<base::OperationMultipleEval> B;

 public:
  /**
   * Std-Constructor
   *
   * @param grid reference to the sparse grid
   * @param lhsData reference to base::DataMatrix that contains the denominator
   * data points
   * @param rhsData reference to base::DataMatrix that contains the numerator
   * data points
   * @param C the regression functional
   * @param lambdaRegression the lambda, the regression parameter
   */
  DMSystemMatrixDRE(base::Grid& grid, base::DataMatrix& lhsData,
                    sgpp::base::DataMatrix& rhsData,
                    std::shared_ptr<base::OperationMatrix> C,
                    double lambdaRegression);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixDRE();

  virtual void mult(base::DataVector& alpha, base::DataVector& result);

  /**
   * Generates the right hand side of the classification equation
   *
   * @param classes the class information of the training data
   * @param b reference to the vector that will contain the result of the matrix
   * vector
   *   multiplication on the rhs
   */
  virtual void generateb(base::DataVector& b);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DMSYSTEMMATRIXDRE_HPP */
