// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

class DBMatOnlineDE_SMW : public DBMatOnlineDE {
 public:
  /**
   * Constructor
   *
   * @param offline The offline object we base our evaluations on.
   * @param lambda The regularization strength
   * @param grid The underlying grid
   * @param beta The initial weighting factor
   */
  explicit DBMatOnlineDE_SMW(DBMatOffline& offline, Grid& grid, double lambda,
                             double beta = 0.);

  /**
   * Returns the additive component of the sherman-morrison-formula, which
   * yields all the information about the refined points
   */
  sgpp::base::DataMatrix& getB() { return this->b_adapt_matrix_; }

  /**
   * @returns distributed version of matrix B
   */
  DataMatrixDistributed& getBDistributed() {
    return this->b_adapt_matrix_distributed_;
  }

  /**
   * Adds new DataVector to list of refined points
   * For testing purposes only
   *
   * @param x The DataVector to add
   */
  void add_new_refine_point(sgpp::base::DataVector& x) {
    this->refined_points_.push_back(x);
  }

  /**
   * Gets pointer to the container of refined points, only for testing purposes
   */
  std::vector<sgpp::base::DataVector>* getRefinedPointsPointer() {
    return &(this->refined_points_);
  }

  /**
   * Rank-one updates/downdates the system matrix, based on the
   * Sherman-Morrison-formula
   * In the current version of the function, the refinePts already are adapted
   * to the regularization parameter lambda.
   *
   * @param X Refine/Coarsen Matrix of L2-Products
   * @param newPoints number of refined points
   * @param refine decides: true for refining, false for coarsening
   * @param coarsen_indices the indices of points to coarsen
   */
  void smw_adapt(DataMatrix& X, size_t newPoints, bool refine,
                 std::vector<size_t> coarsen_indices = {});

  /**
   * Parallel/Distributed version of smw_adapt, uses Scalapack
   *
   * @param X Refine/Coarsen Matrix of L2-Products
   * @param newPoints number of refined points
   * @param refine decides: true for refining, false for coarsening
   * @param coarsen_indices the indices of points to coarsen
   * @param processGrid
   * @param parallelConfig
   */
  void smw_adapt_parallel(DataMatrixDistributed& X, size_t newPoints,
                          bool refine,
                          std::shared_ptr<BlacsProcessGrid> processGrid,
                          const ParallelConfiguration& parallelConfig,
                          std::vector<size_t> coarsen_indices = {});

  /**
   * @param densityEstimationConfig configuration for the density estimation
   * @param grid the underlying grid
   * @param numAddedGridPoints Number of grid points inserted at the end of the
   * grid storage
   * @param deletedGridPointIndices Indices of grid points that were deleted
   * @param lambda The last best lambda value
   * @return list of grid points, that cannot be coarsened
   */
  std::vector<size_t> updateSystemMatrixDecomposition(
      DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
      size_t numAddedGridPoints, std::vector<size_t>& deletedGridPointIndices,
      double lambda) override;

  /**
   * Parallel/Distributed version of updateSystemMatrixDecomposition, uses
   * scalapack and
   * needs to be called with additional parameters for distributed matrices
   * @param densityEstimationConfig configuration for the density estimation
   * @param grid the underlying grid
   * @param numAddedGridPoints Number of grid points inserted at the end of the
   * grid storage
   * @param deletedGridPointIndices Indices of grid points that were deleted
   * @param lambda The last best lambda value
   * @param processGrid
   * @param parallelConfig
   * @return list of grid points, that cannot be coarsened
   */
  std::vector<size_t> updateSystemMatrixDecompositionParallel(
      DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
      size_t numAddedGridPoints, std::vector<size_t>& deletedGridPointIndices,
      double lambda, std::shared_ptr<BlacsProcessGrid> processGrid,
      const ParallelConfiguration& parallelConfig);

  /**
   * Synchronizes the distributed decomposition, only has an effect if scalapack
   * is used.
   */
  void syncDistributedDecomposition(
      std::shared_ptr<BlacsProcessGrid> processGrid,
      const ParallelConfiguration& parallelConfig) override;

 protected:
  // matrix, which holds information about refined/coarsened points
  sgpp::base::DataMatrix b_adapt_matrix_;

  // distributed version of matrix B, only initialized if scalapack is used
  DataMatrixDistributed b_adapt_matrix_distributed_;

  // holds all prior refined points, to know what to coarsen later on
  std::vector<sgpp::base::DataVector> refined_points_;

  // points to end of refined_points_, !which already were processed!
  size_t current_refine_index;

  // tells if refinement/coarsening has been performed on the matrix
  // b_adapt_matrix_ yet
  bool b_is_refined;

  /**
   * Solves the system (R + lambda*I) * alpha = b, and obtains alpha
   * The solving is done after offline and online phase and works as follows:
   * (R + lambda*I)^{-1} * b = ((Q*T^{-1}*Q^{t} + B) * b = alpha,
   * where B holds the refine/coarsen information
   *
   * @param alpha the datavetor where surplusses will be stored
   * @param b The right hand side of the system
   * @param grid the underlying grid
   * @param densityEstimationConfig configuration for the density estimation
   * @param do_cv Specifies, if cross-validation should be done (todo: currently
   * not implemented)
   */
  void solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
                DensityEstimationConfiguration& densityEstimationConfig,
                bool do_cv) override;

  /**
   * Solves the distributed system (R + lambda*I) * alpha = b in parallel and
   * obtains alpha.
   * For more information, see solveSLE.
   *
   * @param alpha distributed datavector for surplusses
   * @param b distributed right hand side
   * @param grid the underlying grid
   * @param densityEstimationConfig configuration for the density estimation
   * @param do_cv cross-validation (currently not implemented)
   */
  void solveSLEParallel(DataVectorDistributed& alpha, DataVectorDistributed& b,
                        Grid& grid,
                        DensityEstimationConfiguration& densityEstimationConfig,
                        bool do_cv) override;

 private:
  /**
   * Computes the L_2 products of the refined gridpoints and pushes them into
   * the
   * refined_points_ container member.
   *
   * @param X the matrix for storing the computed values
   * @param grid the underlying grid
   * @param newPoints the number of points to refine
   * @param newLambda the regularization coefficient added to the diagonal
   * elements
   */
  void compute_L2_refine_matrix(DataMatrix& X, Grid& grid, size_t newPoints,
                                double newLambda);

  /**
   * Computes the L_2 products of the gridpoints to coarsen.
   * Note: member refined_points_ doesn't get changed here, but in the
   * corresponding smw_adapt call
   *
   * @param X the matrix for storing the computed values
   * @param grid the underlying grid
   * @param coarsen_indices the indices of points to coarsen
   */
  void compute_L2_coarsen_matrix(DataMatrix& X, Grid& grid,
                                 std::vector<size_t> coarsen_indices);
};
}  // namespace datadriven
}  // namespace sgpp
