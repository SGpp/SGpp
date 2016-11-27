/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * IChol.cpp
 *
 *  Created on: Nov 26, 2016
 *      Author: Michael Lettrich
 */

#include "IChol.hpp"

namespace sgpp {
namespace datadriven {

void IChol::decompose(DataMatrix& matrix, size_t sweeps) {
  // for all sweeps
  for (auto sweep = 0u; sweep < sweeps; sweep++) {
    // for each row
    for (auto i = 0u; i < matrix.getNrows(); i++) {
      // in each column until diagonal element
      for (auto j = 0u; j <= i; j++) {
        // calculate sum;
        auto s = matrix.get(i, j);
        for (auto k = 0u; k < j; k++) {
          s -= matrix.get(i, k) * matrix.get(j, k);
        }

        // Update value;
        if (i != j) {
          matrix.set(i, j, s / matrix.get(j, j));
        } else {
          matrix.set(j, j, sqrt(s));
        }
      }
    }
  }
}

void IChol::decompose(const DataMatrix& matrix, DataMatrix& result, size_t sweeps) {
  result = matrix;
  decompose(result, sweeps);
}

void IChol::updateLastNRows(DataMatrix& matrix, size_t numRows, size_t sweeps) {}

} /* namespace datadriven */
} /* namespace sgpp */
