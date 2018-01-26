/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineIGE.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization as a base class for all algorithms based on gaussian elimination on
 * a dense matrix.
 */
class DBMatOfflineGE : public DBMatOffline {
 public:
  explicit DBMatOfflineGE(const sgpp::base::RegularGridConfiguration& gridConfig,
                          const sgpp::base::AdpativityConfiguration& adaptivityConfig,
                          const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                          const sgpp::datadriven::DecompositionConfiguration& decompositionConfig);

  explicit DBMatOfflineGE(const std::string& fileName);

  /**
   * Builds the right hand side matrix with identity regularization term
   */
  void buildMatrix() override;

 protected:
  DBMatOfflineGE();
};

} /* namespace datadriven */
} /* namespace sgpp */
