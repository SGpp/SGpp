/*
 * DBMatOfflineLU.hpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOfflineLU : public DBMatOfflineGE {
 public:
  DBMatOfflineLU(DBMatDensityConfiguration& oc);

  void decomposeMatrix() override;

  void permuteVector(DataVector& b);
};

} /* namespace datadriven */
} /* namespace sgpp */
