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

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

/**
 * Class that implements the virtual class SGPP::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
class DMSystemMatrixVectorizedIdentity : public
  SGPP::datadriven::DMSystemMatrixBase {
 private:
  /// vectorization mode
  VectorizationType vecMode_;
  /// Number of orignal training instances
  size_t numTrainingInstances_;
  /// Number of patched and used training instances
  size_t numPatchedTrainingInstances_;
  /// OperationB for calculating the data matrix
  SGPP::parallel::OperationMultipleEvalVectorized* B_;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to SGPP::base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   * @param vecMode vectorization mode
   */
  DMSystemMatrixVectorizedIdentity(SGPP::base::Grid& SparseGrid,
                                   SGPP::base::DataMatrix& trainData, double lambda, VectorizationType vecMode);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixVectorizedIdentity();

  virtual void mult(SGPP::base::DataVector& alpha,
                    SGPP::base::DataVector& result);

  virtual void generateb(SGPP::base::DataVector& classes,
                         SGPP::base::DataVector& b);

  virtual void rebuildLevelAndIndex();
};

}
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP */