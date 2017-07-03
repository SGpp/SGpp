/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEFactory.hpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

namespace sgpp {
namespace datadriven {
namespace DBMatOnlineDEFactory {

/**
 * Factory to build a DBMatOnlineDE object to manipulate the decomposition in offline object
 * @param offline offline object that holds the decomposed system matrix
 * @param beta plasticity weighting factor. If set to 0, no plasticity takes place.
 */
DBMatOnlineDE* buildDBMatOnlineDE(DBMatOffline& offline, double beta = 0);
} /* namespace DBMatOnlineDEFactory */
} /* namespace datadriven */
} /* namespace sgpp */
