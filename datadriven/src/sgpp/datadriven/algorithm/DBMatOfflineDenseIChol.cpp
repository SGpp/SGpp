/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineDenseIChol.cpp
 *
 *  Created on: Apr 15, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineDenseIChol::DBMatOfflineDenseIChol(const DBMatDensityConfiguration& oc)
    : DBMatOfflineChol(oc) {}

DBMatOfflineDenseIChol::DBMatOfflineDenseIChol(const std::string& fileName)
    : DBMatOfflineChol{fileName} {}

DBMatOffline* DBMatOfflineDenseIChol::clone() { return new DBMatOfflineDenseIChol{*this}; }

void DBMatOfflineDenseIChol::decomposeMatrix() {
  if (isConstructed) {
    if (isDecomposed) {
      return;
    } else {
      // extract lower triangular matrix.
      for (size_t i = 0; i < lhsMatrix.getNrows() - 1; i++) {
        for (size_t j = i + 1; j < lhsMatrix.getNcols(); j++) {
          lhsMatrix.set(i, j, 0);
        }
      }
      IChol::decompose(lhsMatrix, 1);
    }
    isDecomposed = true;
  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
