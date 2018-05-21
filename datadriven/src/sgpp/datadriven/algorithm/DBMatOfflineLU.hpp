/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineLU.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */

#ifdef USE_GSL
#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

#include <gsl/gsl_permutation.h>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization that uses a LU factorization on
 * a dense matrix. Does not allow refinement nor changes of regularization parameter.
 */
class DBMatOfflineLU : public DBMatOfflineGE {
 public:
  explicit DBMatOfflineLU(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::base::AdpativityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

  explicit DBMatOfflineLU(const std::string& fileName);

  DBMatOfflineLU(const DBMatOfflineLU& rhs);

  DBMatOfflineLU(DBMatOfflineLU&& rhs) = default;

  DBMatOfflineLU& operator=(const DBMatOfflineLU& rhs);

  DBMatOfflineLU& operator=(DBMatOfflineLU&& rhs) = default;

  virtual ~DBMatOfflineLU() = default;

  DBMatOffline* clone() override;

  /**
   * This decomposition type is not refineable.
   * @return always returns false;
   */
  bool isRefineable() override;

  void decomposeMatrix() override;

  /**
   * Apply permutation vector to the LU factors
   * @param b permutation vector
   */
  void permuteVector(DataVector& b);

  void store(const std::string& fname) override;

 private:
  /**
   * Stores the permutation that was applied on the matrix during decomposition for stability
   * reasons.
   */
  std::unique_ptr<gsl_permutation> permutation;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* USE_GSL */
