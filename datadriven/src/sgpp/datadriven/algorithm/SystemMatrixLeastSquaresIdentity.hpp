// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSUBSPACES_HPP
#define DMSYSTEMMATRIXSUBSPACES_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

// #include <AbstractOperationMultipleEval.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
class SystemMatrixLeastSquaresIdentity : public datadriven::DMSystemMatrixBase {
 private:
  /// vectorization mode
  // ComputeKernelType kernelType;
  /// Number of orignal training instances
  size_t instances;
  /// Number of patched and used training instances
  size_t paddedInstances;
  /// OperationB for calculating the data matrix
  // AbstractOperationMultipleEval* B;
  std::unique_ptr<base::OperationMultipleEval> B;

  base::Grid& grid;

  datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   */
  SystemMatrixLeastSquaresIdentity(base::Grid& SparseGrid, base::DataMatrix& trainData,
                                   double lambda);

  /**
   * Std-Destructor
   */
  virtual ~SystemMatrixLeastSquaresIdentity();

  virtual void mult(base::DataVector& alpha, base::DataVector& result);

  virtual void generateb(base::DataVector& classes, base::DataVector& b);

  virtual void prepareGrid();

  void setImplementation(datadriven::OperationMultipleEvalConfiguration operationConfiguration);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DMSYSTEMMATRIXSUBSPACES_HPP */
