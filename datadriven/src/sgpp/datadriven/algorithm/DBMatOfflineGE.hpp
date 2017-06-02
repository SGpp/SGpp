/*
 * DBMatOfflineIGE.hpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization as a base class for all algorithms based on gaussian elimination on
 * a dense matrix.
 */
class DBMatOfflineGE : public DBMatOffline {
 public:
  explicit DBMatOfflineGE(const DBMatDensityConfiguration& oc);

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
