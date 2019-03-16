/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEChol.hpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

#include <list>
#include <vector>

/**
 * Class that stores, generates and manipulates a density function during online phase in on/off
 * learning. This specialization operates on offline objects based on different Cholesky
 * decompositions.
 */
namespace sgpp {
namespace datadriven {

using sgpp::datadriven::DataVector;

class DBMatOnlineDEChol : public DBMatOnlineDE {
 public:
  /**
   * Constructor
   *
   * @param offline The offline object we base our evaluations on.
   * @param lambda The regularization strength (TODO(fuchsgruber) remove this)
   * @param grid The underlying grid (TODO(fuchsgruber) do we need this?)
   * @param beta The initial weighting factor
   */
  explicit DBMatOnlineDEChol(DBMatOffline& offline, Grid& grid, double lambda, double beta = 0.);

  /**
   * Delegates call to choleskyModification
   * @param densityEstimationConfig Configuration to the density estimation
   * @param grid the underlying grid
   * @param numAddedGridPoints Number of grid points inserted at the end of the grid storage
   * @param deletedGridPointIndices Indices of grid points that were deleted
   * @param lambda The last best lambda value
   * @return list of grid points, that cannot be coarsened
   */
  std::vector<size_t> updateSystemMatrixDecomposition(
      DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
      size_t numAddedGridPoints, std::list<size_t> deletedGridPointIndices, double lambda) override;

 protected:
  void solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
                DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) override;

  /**
   * Not implemented for this decomposition
   */
  void solveSLEParallel(DataVectorDistributed& alpha, DataVectorDistributed& b, Grid& grid,
                        DensityEstimationConfiguration& densityEstimationConfig,
                        const ParallelConfiguration& parallelConfig,
                        std::shared_ptr<BlacsProcessGrid> processGrid, bool do_cv) override {
    // TODO(jan) implement for this decomposition
    throw base::not_implemented_exception(
        "Distributed parallel solve not implemented for this decomposition");
  }

  DBMatDMSChol* buildCholSolver(DBMatOffline& offlineObject, Grid& grid,
                                DensityEstimationConfiguration& densityEstimationConfig,
                                bool doCV) const;
};

} /* namespace datadriven */
} /* namespace sgpp */
