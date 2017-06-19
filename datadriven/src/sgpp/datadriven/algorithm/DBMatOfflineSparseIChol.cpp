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

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineSparseIChol.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineSparseIChol::DBMatOfflineSparseIChol(const DBMatDensityConfiguration& oc)
    : DBMatOfflineDenseIChol(oc) {}

DBMatOfflineSparseIChol::DBMatOfflineSparseIChol(const std::string& fileName)
    : DBMatOfflineDenseIChol{fileName} {}

DBMatOffline* sgpp::datadriven::DBMatOfflineSparseIChol::clone() {
  return new DBMatOfflineSparseIChol{*this};
}

void DBMatOfflineSparseIChol::decomposeMatrix() {
  if (isConstructed) {
    if (isDecomposed) {
      return;
    } else {
      //  auto begin = std::chrono::high_resolution_clock::now();

      DataMatrix tmpMatrix{lhsMatrix.getNrows(), lhsMatrix.getNcols()};

// only copy lower triangular matrix
#pragma omp parallel for schedule(guided)
      for (auto i = 0u; i < tmpMatrix.getNrows(); i++) {
#pragma omp simd
        for (auto j = 0u; j <= i; j++) {
          tmpMatrix.set(i, j, lhsMatrix.get(i, j));
        }
      }

      ichol(tmpMatrix, lhsMatrix, config.icholParameters.sweepsDecompose);
    }
    isDecomposed = true;
    //    auto end = std::chrono::high_resolution_clock::now();
    //    std::cout << "IChol decompostition took"
    //              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<
    //              "ms"
    //              << std::endl;
  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
}

void DBMatOfflineSparseIChol::ichol(const DataMatrix& matrix, DataMatrix& result, size_t sweeps,
                                    size_t startRow) {
  SparseDataMatrix sparseResult(result.getNrows(), result.getNcols());
  SparseDataMatrix::fromDataMatrixTriangular(result, sparseResult);

  const auto matSize = sparseResult.getNrows();
  // get the data vector
  auto& matData = sparseResult.getDataVector();
  // get the rows
  const auto& rowPtrs = sparseResult.getRowPtrVector();
  // get the cols
  const auto& colIndices = sparseResult.getColIndexVector();

  // for all sweeps
  for (auto sweep = 0u; sweep < sweeps; sweep++) {
#pragma omp parallel for
    for (auto dataIter = rowPtrs[startRow]; dataIter < matData.size(); dataIter++) {
      const auto col = colIndices[dataIter];
      const auto row = [&rowPtrs, dataIter]() {
        const auto find = std::lower_bound(std::begin(rowPtrs), std::end(rowPtrs), dataIter);
        return static_cast<size_t>(
            ((*find > dataIter || find == std::end(rowPtrs)) ? (find - 1) : find) -
            rowPtrs.begin());
      }();

      auto s = matrix.get(row, col);

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
  SparseDataMatrix::toDataMatrix(sparseResult, result);
}

} /* namespace datadriven */
} /* namespace sgpp */
