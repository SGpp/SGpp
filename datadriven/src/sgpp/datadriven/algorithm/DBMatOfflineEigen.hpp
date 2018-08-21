/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineEigen.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */
#ifdef USE_GSL

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

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

  DBMatOffline* clone() override;

  /**
   * Returns the decomposition type of the DBMatOffline object
   * @return the type of matrix decomposition
   */
  sgpp::datadriven::MatrixDecompositionType getDecompositionType() override;

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
  void decomposeMatrix(RegularizationConfiguration& regularizationConfig,
      DensityEstimationConfiguration& densityEstimationConfig) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */
