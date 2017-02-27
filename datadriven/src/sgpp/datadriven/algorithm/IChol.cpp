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
  DataMatrix A;
  SparseDataMatrix::toDataMatrix(matrix, A);

  const auto matSize = matrix.getNrows();
  // get the data vector
  auto& matData = matrix.getDataVector();
  // get the rows
  const auto& rowPtrs = matrix.getRowPtrVector();
  // get the cols
  const auto& colIndices = matrix.getColIndexVector();

  // for all sweeps
  for (auto sweep = 0u; sweep < sweeps; sweep++) {
    for (auto dataIter = 0u; dataIter < matData.size(); dataIter++) {
      const auto col = colIndices[dataIter];
      const auto find = std::lower_bound(std::begin(rowPtrs), std::end(rowPtrs), dataIter);
      auto row =
          ((*find > dataIter || find == std::end(rowPtrs)) ? (find - 1) : find) - rowPtrs.begin();

      printf("starting decomp of [%d,%d]= %f \n", row, col, matData[dataIter]);

      auto s = matData[dataIter];

      auto upperFirst = colIndices.begin() + rowPtrs[col];
      const auto upperLast =
          colIndices.begin() + (col + 1 < matSize ? (rowPtrs[col + 1] - 1) : dataIter);
      auto lowerFist = colIndices.begin() + rowPtrs[row];
      const auto lowerLast = colIndices.begin() + dataIter;

      printf("\tBounds:\n\tlower: %d, %d; upper: %d, %d\n", upperFirst - colIndices.begin(),
             upperLast - colIndices.begin(), lowerFist - colIndices.begin(),
             lowerLast - colIndices.begin());

      // sparse dot product by merging in O(n+m)
      while (lowerFist != lowerLast) {
        // if we're out of nonzeors in the upper row, then we're also done.

        if (upperFirst == upperLast) {
          std::cout << "done\n";
          break;
        }
        if (*upperFirst < *lowerFist) {
          // printf("Mismatch: upper %d < %d lower: incrementing to ", *upperFirst, *lowerFist);
          ++upperFirst;
          // std::cout << *upperFirst << std::endl;
        } else if (*lowerFist < *upperFirst) {
          // printf("Mismatch: upper %d < %d lower: incrementing to ", *upperFirst, *lowerFist);
          ++lowerFist;
          // std::cout << *lowerFist << std::endl;
        } else {
          // printf("Match: upper %d == %d lower\n", *upperFirst, *lowerFist);
          printf("\tcalculating [%d,%d]*[%d,%d] => %f-%f*%f=", col, *upperFirst, row, *lowerFist, s,
                 matData[upperFirst - colIndices.begin()], matData[lowerFist - colIndices.begin()]);
          s -= matData[upperFirst - colIndices.begin()] * matData[lowerFist - colIndices.begin()];
          std::cout << s << std::endl;
          ++upperFirst;
          ++lowerFist;
        }
      }

      if (row != col) {
        const auto index = rowPtrs[col + 1] - 1;
        matData[dataIter] = s / matData[index];
        std::cout << "\t" << matData[dataIter] << " = " << s << " / " << matData[index]
                  << std::endl;
      } else {
        matData[dataIter] = sqrt(s);
        std::cout << "\t" << matData[dataIter] << " = "
                  << "sqrt(" << s << ")" << std::endl;
      }
    }

    DataMatrix A;
    SparseDataMatrix::toDataMatrix(matrix, A);
    std::cout << "after sweep " << sweep << ":\n" << A.toString() << "\n";
  }

  // #pragma omp parallel
  //  { /* omp parallel */
  // #pragma omp for

  for (auto i = 0u; i < A.getNrows(); i++) {
    // in each column until diagonal element
    for (auto j = 0u; j < i; j++) {
      if (A.get(i, j) > 0) {
        printf("starting decomp of [%d,%d]= %f \n", i, j, A.get(i, j));
        // calculate sum;
        auto s = A.get(i, j);
        for (auto k = 0u; k < j; k++) {
          if (A.get(i, k) > 0 && A.get(j, k) > 0) {
            printf("\tcalculating [%d,%d]*[%d,%d] => %f-%f*%f=", i, k, j, k, s, A.get(i, k),
                   A.get(j, k));
            s -= A.get(i, k) * A.get(j, k);
            std::cout << s << std::endl;
          }
        }

        A.set(i, j, s / A.get(j, j));
        std::cout << "\t" << A.get(i, j) << " = " << s << " / " << A.get(j, j) << std::endl;
      }
    }

    // do the diagonal element:
    // calculate sum;

    printf("starting decomp of [%d,%d]= %f \n", i, i, A.get(i, i));
    auto s = A.get(i, i);
    for (auto k = 0u; k < i; k++) {
      if (A.get(i, k) > 0) {
        printf("\tcalculating [%d,%d]*[%d,%d] => %f-%f*%f=", i, k, i, k, s, A.get(i, k),
               A.get(i, k));
        s -= A.get(i, k) * A.get(i, k);
        std::cout << s << std::endl;
      }
    }
    A.set(i, i, sqrt(s));
    std::cout << "\t" << A.get(i, i) << " = "
              << "sqrt(" << s << ")" << std::endl;
  }

  std::cout << "Full:\n" << A.toString() << "\n";
  //  } /* omp parallel */
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
