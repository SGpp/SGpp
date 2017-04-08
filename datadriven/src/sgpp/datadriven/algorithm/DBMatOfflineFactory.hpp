/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineFactory.hpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

namespace sgpp {
namespace datadriven {
namespace DBMatOfflineFactory {

DBMatOffline* buildOfflineObject(const DBMatDensityConfiguration& config);

DBMatOffline* buildFromFile(const std::string& fname);

} /* namespace DBMatOfflineFactory */
} /* namespace datadriven */
} /* namespace sgpp */
