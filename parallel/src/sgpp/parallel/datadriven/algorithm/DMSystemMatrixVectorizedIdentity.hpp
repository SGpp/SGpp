// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace parallel {

/**
 * Class that implements the virtual class sgpp::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
class DMSystemMatrixVectorizedIdentity : public sgpp::datadriven::DMSystemMatrixBase {
 private:
  /// vectorization mode
  VectorizationType vecMode_;
  /// Number of orignal training instances
  size_t numTrainingInstances_;
  /// Number of patched and used training instances
  size_t numPatchedTrainingInstances_;
  /// OperationB for calculating the data matrix
  std::unique_ptr<sgpp::parallel::OperationMultipleEvalVectorized> B_;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to sgpp::base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   * @param vecMode vectorization mode
   */
  DMSystemMatrixVectorizedIdentity(sgpp::base::Grid& SparseGrid, sgpp::base::DataMatrix& trainData,
                                   double lambda, VectorizationType vecMode);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixVectorizedIdentity();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void generateb(sgpp::base::DataVector& classes, sgpp::base::DataVector& b);

  virtual void prepareGrid() { rebuildLevelAndIndex(); }

  virtual void rebuildLevelAndIndex();
};

}  // namespace parallel
}  // namespace sgpp

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP */
