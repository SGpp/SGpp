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

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * factories to build the specialization of the DBMatOffline objects.
 */
namespace DBMatOfflineFactory {

/**
 * Based on the configuration file, build the appropriate DBMatOffline object and return it
 * @param config Configuration file with parameters set
 * @return new instance of DBMatOffline implementor owned by caller.
 */
DBMatOffline* buildOfflineObject(const DBMatDensityConfiguration& config);

/**
 * Read a serialized DBMatOffline object and construct a new object with the information.
 * @param fname Path to the serialized DBMatOffline object.
 * @return new instance of DBMatOffline implementor owned by caller.
 */
DBMatOffline* buildFromFile(const std::string& fname);

} /* namespace DBMatOfflineFactory */
} /* namespace datadriven */
} /* namespace sgpp */
