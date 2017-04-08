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

class DBMatOfflineGE : public DBMatOffline {
 public:
  explicit DBMatOfflineGE(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineGE(const std::string& fileName);

  /**
   * Builds the right hand side matrix with or without the regularization term
   * depending
   * on the type of decomposition
   */
  void buildMatrix() override;

 protected:
  DBMatOfflineGE();
};

} /* namespace datadriven */
} /* namespace sgpp */
