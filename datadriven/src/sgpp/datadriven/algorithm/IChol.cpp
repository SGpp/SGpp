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

#include <cmath>
#include <iostream>

namespace sgpp {
namespace datadriven {

void IChol::decompose(SparseDataMatrix& matrix, size_t sweeps) {
  const auto matSize = matrix.getNrows();
  // get the data vector
  auto& matData = matrix.getDataVector();
  // get the rows
  const auto& rowPtrs = matrix.getRowPtrVector();
  // get the cols
  const auto& colIndices = matrix.getColIndexVector();

  // for all sweeps
  for (auto sweep = 0u; sweep < sweeps; sweep++) {
    //#pragma omp parallel for
    for (auto dataIter = 0u; dataIter < matData.size(); dataIter++) {
      const auto col = colIndices[dataIter];
      const auto row = [&]() {
        const auto find = std::lower_bound(std::begin(rowPtrs), std::end(rowPtrs), dataIter);
        return static_cast<size_t>(
            ((*find > dataIter || find == std::end(rowPtrs)) ? (find - 1) : find) -
            rowPtrs.begin());
      }();

      auto s = matData[dataIter];

      auto upperFirst = colIndices.begin() + rowPtrs[col];
      const auto upperLast =
          colIndices.begin() + (col + 1 < matSize ? (rowPtrs[col + 1] - 1) : dataIter);
      auto lowerFist = colIndices.begin() + rowPtrs[row];
      const auto lowerLast = colIndices.begin() + dataIter;

      // sparse dot product by merging in O(n+m)
      while (lowerFist != lowerLast) {
        // if we're out of nonzeors in the upper row, then we're also done.
        if (upperFirst == upperLast) {
          break;
        }
        if (*upperFirst < *lowerFist) {
          ++upperFirst;
        } else if (*lowerFist < *upperFirst) {
          ++lowerFist;
        } else {
          s -= matData[upperFirst - colIndices.begin()] * matData[lowerFist - colIndices.begin()];
          ++upperFirst;
          ++lowerFist;
        }
      }

      if (row != col) {
        const auto index = rowPtrs[col + 1] - 1;
        matData[dataIter] = s / matData[index];
      } else {
        matData[dataIter] = sqrt(s);
      }
    }
  }
}

void IChol::updateLastNRows(SparseDataMatrix& matrix, size_t numRows, size_t sweeps) {}

void IChol::normToUnitDiagonal(SparseDataMatrix& matrix, DataVector& norms) {
  const auto matSize = matrix.getNrows();
  norms.resize(matSize);

  // get the data vector
  auto& matData = matrix.getDataVector();
  // get the rows
  const auto& rowPtrs = matrix.getRowPtrVector();
  // get the cols
  const auto& colIndices = matrix.getColIndexVector();

#pragma omp parallel for
  for (auto i = 0u; i < matSize - 1; i++) {
    // for a tridiagonal matrix, the last element in a row is the diagonal element:
    const auto rowIdx = rowPtrs[i + 1] - 1;
    norms[i] = 1.0 / sqrt(matData[rowIdx]);
  }
  // peel of last element to avoid branching for better auto vectorization.
  norms[matSize - 1] = 1.0 / sqrt(matData.back());

// Calculate (D*A*D)
#pragma omp parallel for
  for (auto i = 0u; i < matSize - 1; i++) {
    const auto leftVecVal = norms[i];
    for (auto j = rowPtrs[i]; j < rowPtrs[i + 1]; j++) {
      const auto rightVecVal = norms[colIndices[j]];
      const auto res = leftVecVal * matData[j] * rightVecVal;
      matData[j] = res;
    }
  }

  const auto i = matSize - 1;
  // peel of last element to avoid branching for better auto vectorization.
  for (auto j = rowPtrs[i]; j < matData.size(); j++) {
    const auto rightVecVal = norms[colIndices[j]];
    const auto res = norms[i] * matData[j] * rightVecVal;
    matData[j] = res;
  }
}

void IChol::reaplyDiagonal(SparseDataMatrix& matrix, DataVector& norms) {
  const auto matSize = matrix.getNrows();
  // get the data vector
  auto& matData = matrix.getDataVector();
  // get the rows
  const auto& rowPtrs = matrix.getRowPtrVector();
  // get the cols
  const auto& colIndices = matrix.getColIndexVector();

#pragma omp parallel for
  for (auto i = 0u; i < matSize - 1; i++) {
    const auto leftVecVal = norms[i];
    for (auto j = rowPtrs[i]; j < rowPtrs[i + 1]; j++) {
      const auto rightVecVal = norms[colIndices[j]];
      const auto res = matData[j] / (leftVecVal * rightVecVal);
      matData[j] = res;
    }
  }

  const auto i = matSize - 1;
  for (auto j = rowPtrs[i]; j < matData.size(); j++) {
    const auto rightVecVal = norms[colIndices[j]];
    const auto res = matData[j] / (norms[i] * rightVecVal);
    matData[j] = res;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
