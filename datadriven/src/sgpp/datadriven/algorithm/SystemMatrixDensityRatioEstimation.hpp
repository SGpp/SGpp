// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixDRE.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

// #include <AbstractOperationMultipleEval.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements the virtual class base::OperationMatrix for the application of density
 * ratio estimation for SystemMatrixDRE
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions vectorized formulations are used.
 */
class SystemMatrixDensityRatioEstimation : public datadriven::DMSystemMatrixDRE {
 private:
  /// vectorization mode
  // ComputeKernelType kernelType;
  /// Number of original training instances of first dataset
  size_t instancesP;
  /// Number of patched and used training instances of first dataset
  size_t paddedInstancesP;
  /// Number of original training instances of second dataset
  size_t instancesQ;
  /// Number of patched and used training instances of second dataset
  size_t paddedInstancesQ;
  /// OperationB for calculating the data matrix for the two datasets
  // AbstractOperationMultipleEval* B;
  std::unique_ptr<base::OperationMultipleEval> B_p;
  std::unique_ptr<base::OperationMultipleEval> B_q;

  base::Grid& grid;

  datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainDataP reference to base::DataMatrix that contains the training data for first
   * dataset
   * @param trainDataQ reference to base::DataMatrix that contains the training data for second
   * dataset
   * @param lambda the lambda, the regression parameter
   */
  SystemMatrixDensityRatioEstimation(base::Grid& SparseGrid, base::DataMatrix& trainDataP,
                                     base::DataMatrix& trainDataQ, double lambda);

  /**
   * Std-Destructor
   */
  virtual ~SystemMatrixDensityRatioEstimation();

  virtual void mult(base::DataVector& alpha, base::DataVector& result);

  virtual void generateb(base::DataVector& b);

  virtual void prepareGrid();

  void setImplementation(datadriven::OperationMultipleEvalConfiguration operationConfiguration);
};

}  // namespace datadriven
}  // namespace sgpp
