/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineFactory.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineIChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>

namespace sgpp {
namespace datadriven {
using sgpp::base::algorithm_exception;

DBMatOffline* DBMatOfflineFactory::buildOfflineObject(const DBMatDensityConfiguration& config) {
  auto type = config.decomp_type_;

  switch (type) {
    case (DBMatDecompostionType::DBMatDecompEigen):
      return new DBMatOfflineEigen(config);
      break;
    case (DBMatDecompostionType::DBMatDecompLU):
      return new DBMatOfflineLU(config);
      break;
    case (DBMatDecompostionType::DBMatDecompChol):
      return new DBMatOfflineChol(config);
      break;
    case (DBMatDecompostionType::DBMatDecompIChol):
      return new DBMatOfflineIChol(config);
      break;
    default:
      throw algorithm_exception("Trying to build offline object from unknown decomposition type");
      return nullptr;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
