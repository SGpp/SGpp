// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

/**
 * Class that stores, generates and manipulates a density function during online phase in on/off
 * learning. This specialization operates on offline objects based on a LU decomposition.
 */
class DBMatOnlineDELU : public DBMatOnlineDE {
 public:
  /**
   * Constructor
   *
   * @param offline The offline object we base our evaluations on.
   * @param lambda The regularization strength (TODO(fuchsgruber) remove this)
   * @param grid The underlying grid (TODO(fuchsgruber) do we need this?)
   * @param beta The initial weighting factor
   */
  explicit DBMatOnlineDELU(DBMatOffline& offline, Grid& grid, double lambda, double beta = 0.);

 protected:
  /**
   * Solves the SLE for this matrix decomposition
   * @param alpha the surplusses
   * @param b the rhs of the system
   * @param grid the underlying grid
   * @param densityEstimationConfig configuration for the density estimation
   * @param do_cv whether cross validation should be performed
   */
  void solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
                DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) override;

  /**
   * Not implemented for this decomposition
   */
  void solveSLEParallel(DataVectorDistributed& alpha, DataVectorDistributed& b, Grid& grid,
                        DensityEstimationConfiguration& densityEstimationConfig,
                        bool do_cv) override {
    throw base::not_implemented_exception(
        "Distributed parallel solve not implemented for this decomposition");
  }
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /*USE_GSL*/
