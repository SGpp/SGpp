// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP

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
class DMSystemMatrixVectorizedIdentityMPI : public sgpp::datadriven::DMSystemMatrixBase {
 private:
  /// vectorization mode
  VectorizationType vecMode_;
  /// Number of orignal training instances
  size_t numTrainingInstances_;
  /// Number of patched and used training instances
  size_t numPatchedTrainingInstances_;
  /// OperationB for calculating the data matrix
  std::unique_ptr<sgpp::parallel::OperationMultipleEvalVectorized> B_;

  double waitting_time;

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param trainData reference to sgpp::base::DataMatrix that contains the training data
   * @param lambda the lambda, the regression parameter
   * @param vecMode vectorization mode
   */
  DMSystemMatrixVectorizedIdentityMPI(sgpp::base::Grid& SparseGrid,
                                      sgpp::base::DataMatrix& trainData, double lambda,
                                      VectorizationType vecMode);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixVectorizedIdentityMPI();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void generateb(sgpp::base::DataVector& classes, sgpp::base::DataVector& b);

  virtual void rebuildLevelAndIndex();

  virtual void prepareGrid() { rebuildLevelAndIndex(); }

 private:
  /// how to distribute storage array
  int* _mpi_grid_sizes;
  int* _mpi_grid_offsets;

  /// reference to grid. needed to get new grid size after it changes
  sgpp::base::Grid& m_grid;

  /// how to distribute dataset
  int* _mpi_data_sizes;
  int* _mpi_data_offsets;

  /**
   * Wrapper function that handles communication after calculation and time measurement
   */
  void multVec(base::DataVector& alpha, base::DataVector& result);

  /**
   * Wrapper function that handles communication after calculation and time measurement
   */
  void multTransposeVec(base::DataVector& source, base::DataVector& result);
};

}  // namespace parallel
}  // namespace sgpp

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITYMPI_HPP */
