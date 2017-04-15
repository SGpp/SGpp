/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineIChol.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineIChol.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineIChol::DBMatOfflineIChol(const DBMatDensityConfiguration& oc) : DBMatOfflineChol(oc) {}

DBMatOfflineIChol::DBMatOfflineIChol(const std::string& fileName) : DBMatOfflineChol{fileName} {}

DBMatOffline* sgpp::datadriven::DBMatOfflineIChol::clone() { return new DBMatOfflineIChol{*this}; }

void DBMatOfflineIChol::decomposeMatrix() {
  //  std::cout << "ichol decomposing\n";
  if (isConstructed) {
    if (isDecomposed) {
      return;
    } else {
      //      std::cout << "Full:\n" << lhsMatrix.toString() << "\n\n";
      // extract lower triangular matrix.
      for (size_t i = 0; i < lhsMatrix.getNrows() - 1; i++) {
        for (size_t j = i + 1; j < lhsMatrix.getNcols(); j++) {
          lhsMatrix.set(i, j, 0);
        }
      }

      //      DataMatrix tmp;
      SparseDataMatrix sparseLHS(lhsMatrix);
      //      std::cout << "DataMatrix:\n" << lhsMatrix.toString() << "\n\n";

      DataVector norm{lhsMatrix.getNrows()};
      IChol::normToUnitDiagonal(sparseLHS, norm);
      //      SparseDataMatrix::toDataMatrix(sparseLHS, tmp);
      //      std::cout << "Normed:\n" << tmp.toString() << "\n\n";

      IChol::decompose(sparseLHS, 1);
      //      SparseDataMatrix::toDataMatrix(sparseLHS, tmp);
      //      std::cout << "decomposed:\n" << tmp.toString() << "\n\n";

      IChol::reaplyDiagonal(sparseLHS, norm);
      SparseDataMatrix::toDataMatrix(sparseLHS, lhsMatrix);
      //      std::cout << "reapplied:\n" << lhsMatrix.toString() << "\n\n";
    }

    isDecomposed = true;

  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
