// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

#include <list>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

/**
 * DBMatOffline specialization that uses a cholesky factorization on
 * a dense matrix. The resulting factorization can be updated when the grid changes.
 */
class DBMatOfflineChol : public DBMatOfflineGE {
 public:
  DBMatOfflineChol();

  explicit DBMatOfflineChol(const std::string& fileName);

  DBMatOffline* clone() const override;

  bool isRefineable() override;

  /**
   * Returns the decomposition type of the DBMatOffline object
   * @return the type of matrix decomposition
   */
  sgpp::datadriven::MatrixDecompositionType getDecompositionType() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   * @param regularizationConfig the regularization configuration
   * @param densityEstimationConfig the density estimation configuration
   */
  void decomposeMatrix(const RegularizationConfiguration& regularizationConfig,
                       const DensityEstimationConfiguration& densityEstimationConfig) override;

  /**
   * The parallel/distributed version of decomposeMatrix(...)
   * @param regularizationConfig the regularization configuration
   * @param densityEstimationConfig the density estimation configuration
   * @param processGrid process grid to distribute the matrix on
   * @param parallelConfig
   */
  void decomposeMatrixParallel(RegularizationConfiguration& regularizationConfig,
                               DensityEstimationConfiguration& densityEstimationConfig,
                               std::shared_ptr<BlacsProcessGrid> processGrid,
                               const ParallelConfiguration& parallelConfig) override;

  /**
   * Updates offline cholesky factorization based on coarsed (deletedPoints)
   * and refined (newPoints) gridPoints
   *
   * @param grid the underlying grid
   * @param densityEstimationConfig configuration for the density estimation
   * @param newPoints amount of refined points
   * @param deletedPoints list of indices of last coarsed points
   * @param lambda the regularization parameter
   */
  virtual void choleskyModification(
      Grid& grid, datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      size_t newPoints, std::list<size_t> deletedPoints, double lambda, double c = 0.0);

  /*
   * explicitly computes the inverse
   * note: the computed inverse is not the inverse of the decomposed matrix,
   * but rather the inverse of the previous undecomposed matrix
   * @param inv the matrix to store the computed inverse
   */
  void compute_inverse() override;

  /**
   * parallel/distributed version of compute_inverse()
   * @param processGrid process grid to distribute the matrix on
   * @param parallelConfig
   */
  void compute_inverse_parallel(std::shared_ptr<BlacsProcessGrid> processGrid,
                                const ParallelConfiguration& parallelConfig) override;

 protected:
  /**
   * Permutes the rows of the cholesky factor based on permutations
   * of the system matrix (e.g. coarsening)
   *
   * @param k "left" column to permutate
   * @param l "right" column to permutate
   * @param job = 2        => left circular shift
   *	  1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ...,
   *l-1,l,k,l+1,..., size
   * 	  job = 1       => right circular shift
   * 	  1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ...,
   *l-1,l+1,...size
   */
  void choleskyPermutation(size_t k, size_t l, size_t job);

  /**
   * Updates the cholesky factor when a new grid point is added (e.g. refine)
   *
   * @param newCol DataVector with column to add to the system matrix
   * @param size columns/rows of current Cholesky factor, necessary since the
            allocated memory is increased before the Cholesky factor is modified
   */
  void choleskyAddPoint(DataVector& newCol, size_t size);
};

} /* namespace datadriven */
} /* namespace sgpp */
