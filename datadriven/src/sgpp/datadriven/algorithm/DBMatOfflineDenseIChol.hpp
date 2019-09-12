// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>

#include <list>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

/**
 * DBMatOfflineChol specialization that uses a parallel, iterative incomplete cholesky factorization
 * on a dense matrix. The current implementation is a proof of concept.
 */
class DBMatOfflineDenseIChol : public DBMatOfflineChol {
 public:
  DBMatOfflineDenseIChol();

  explicit DBMatOfflineDenseIChol(const std::string& fileName);

  DBMatOffline* clone() override;

  /**
   * Returns the decomposition type of the DBMatOffline object
   * @return the type of matrix decomposition
   */
  sgpp::datadriven::MatrixDecompositionType getDecompositionType() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   *
   * @param regularizationConfig the regularization configuration
   * @param densityEstimationConfig the density estimation configuration
   */
  void decomposeMatrix(RegularizationConfiguration& regularizationConfig,
      DensityEstimationConfiguration& densityEstimationConfig) override;

  /**
   * Updates offline cholesky factorization based on coarsed (deletedPoints)
   * and refined (newPoints) gridPoints. We ignore coarsening.
   * @param grid the underlying grid
   * @param densityEstimationConfig configuration for the density estimation
   * @param newPoints amount of refined points
   * @param deletedPoints list of indices of last coarsed points that are ignored.
   * @param lambda the regularization parameter
   */
  void choleskyModification(Grid& grid,
      datadriven::DensityEstimationConfiguration& densityEstimationConfig, size_t newPoints,
      std::list<size_t> deletedPoints, double lambda) override;

  /**
   * perform parlallel incomplete cholesky factorization of a matrix. This is an out of place
   * operation.
   * @param matrix the matrix to be decomposed
   * @param result data matrix that will hold the decomposed matrix
   * @param sweeps how many iterations of the algorithm are required until the result is good
   * enough?
   * @param startRow on which row to start the decomposition (needed for refinement)
   */
  static void ichol(const DataMatrix& matrix, DataMatrix& result, size_t sweeps = 4,
                    size_t startRow = 0);
};

} /* namespace datadriven */
} /* namespace sgpp */
