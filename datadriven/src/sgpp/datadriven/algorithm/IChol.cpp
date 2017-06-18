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

void IChol::normToUnitDiagonal(DataMatrix& matrix, DataVector& norms) {
  const auto matSize = matrix.getNrows();
  norms.resize(matSize);
#pragma omp parallel for
  for (auto i = 0u; i < matSize; i++) {
    norms[i] = 1.0 / std::sqrt(matrix.get(i, i));
  }
  std::cout << norms.toString() << std::endl;
// Calculate (D*A*D)
#pragma omp parallel for
  for (auto i = 0u; i < matSize; i++) {
    const auto leftVecVal = norms[i];
    for (auto j = 0u; j <= i; j++) {
      const auto rightVecVal = norms[j];
      const auto res = leftVecVal * matrix.get(i, j) * rightVecVal;
      matrix.set(i, j, res);
    }
  }
}

void IChol::reaplyDiagonal(DataMatrix& matrix, const DataVector& norms) {
#pragma omp parallel for
  for (auto i = 0u; i < matrix.getNrows(); i++) {
    const auto leftVecVal = norms[i];
    for (auto j = 0u; j <= i; j++) {
      const auto rightVecVal = norms[j];
      const auto res = matrix.get(i, j) / (leftVecVal * rightVecVal);
      matrix.set(i, j, res);
    }
  }
}

void IChol::reaplyDiagonalLowerTriangular(DataMatrix& matrix, const DataVector& norms) {
#pragma omp parallel for
  for (auto i = 0u; i < matrix.getNrows(); i++) {
    for (auto j = 0u; j <= i; j++) {
      const auto res = matrix.get(i, j) / norms[i];
      matrix.set(i, j, res);
    }
  }
}

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

void IChol::reaplyDiagonal(SparseDataMatrix& matrix, const DataVector& norms) {
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
  // peel of last element to avoid branching for better auto vectorization.
  for (auto j = rowPtrs[i]; j < matData.size(); j++) {
    const auto rightVecVal = norms[colIndices[j]];
    const auto res = matData[j] / (norms[i] * rightVecVal);
    matData[j] = res;
  }
}

void IChol::reaplyDiagonalLowerTriangular(SparseDataMatrix& matrix, const DataVector& norms) {
  const auto matSize = matrix.getNrows();
  // get the data vector
  auto& matData = matrix.getDataVector();
  // get the rows
  const auto& rowPtrs = matrix.getRowPtrVector();

  //#pragma omp parallel for
  for (auto i = 0u; i < matSize - 1; i++) {
    for (auto j = rowPtrs[i]; j < rowPtrs[i + 1]; j++) {
      matData[j] /= norms[i];
    }
  }

  const auto i = matSize - 1;
  // peel of last element to avoid branching for better auto vectorization.
  for (auto j = rowPtrs[i]; j < matData.size(); j++) {
    matData[j] /= norms[i];
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
