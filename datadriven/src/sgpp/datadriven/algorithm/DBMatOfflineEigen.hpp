// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <gsl/gsl_permutation.h>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization that uses an Eigen factorization on a dense matrix. Eigen
 * factorization allows cheap changing of the regularization parameter but does not allow
 * refinement.
 */
class DBMatOfflineEigen : public DBMatOffline {
 public:
  DBMatOfflineEigen();

  explicit DBMatOfflineEigen(const std::string& fileName);

  DBMatOffline* clone() const override;

  /**
   * Returns the decomposition type of the DBMatOffline object
   * @return the type of matrix decomposition
   */
  sgpp::datadriven::MatrixDecompositionType getDecompositionType() override;

  /**
   * Get the unmodified (without added lambda) system matrix R.
   *
   * @return Matrix R
   */
  const DataMatrix& getUnmodifiedR() override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineEigen::getUnmodifiedR() is not implemented!");
  }

  /**
   * Get the unmodified (without added lambda) system matrix R.
   *
   * @return Matrix R
   */
  const DataMatrixDistributed& getUnmodifiedRDistributed(
      std::shared_ptr<BlacsProcessGrid> processGrid,
      const ParallelConfiguration& parallelConfig) override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineEigen::getUnmodifiedRDistributed() is not implemented!");
  }

  /**
   * Modifies the decomposition to update the regularization parameter lambda
   *
   * @param lambda New lambda value
   */
  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "DBMatOfflineEigen::updateRegularization() is not implemented!");
  }

  /**
   * This decomposition type is not refineable.
   * @return always returns false;
   */
  bool isRefineable() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   *
   * @param regularizationConfig the regularization configuration
   * @param densityEstimationConfig the density estimation configuration
   */
  void decomposeMatrix(const RegularizationConfiguration& regularizationConfig,
                       const DensityEstimationConfiguration& densityEstimationConfig) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */
