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
 * @param grid the underlying grid
 * @param lambda regularization strength (todo(fuchsgruber): maybe remove this)
 * @param beta plasticity weighting factor. If set to 0, no plasticity takes place.
 * @param matDecompType the matrix decomposition type of the online matrix
 */
DBMatOnlineDE* buildDBMatOnlineDE(
    DBMatOffline& offline, Grid& grid, double lambda, double beta = 0,
    MatrixDecompositionType matDecompType = MatrixDecompositionType::Chol);
} /* namespace DBMatOnlineDEFactory */
} /* namespace datadriven */
} /* namespace sgpp */
