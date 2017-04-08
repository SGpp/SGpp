/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDELU.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

sgpp::datadriven::DBMatOnlineDELU::DBMatOnlineDELU(DBMatOffline& offline, double beta)
    : DBMatOnlineDE{offline, beta} {}

void sgpp::datadriven::DBMatOnlineDELU::solveSLE(DataVector& b, bool do_cv) {
  DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

  // Solve the system:
  alpha = DataVector(lhsMatrix.getNcols());
  DBMatDMSBackSub lusolver;
  lusolver.solve(lhsMatrix, alpha, b);
}

} /* namespace datadriven */
} /* namespace sgpp */
