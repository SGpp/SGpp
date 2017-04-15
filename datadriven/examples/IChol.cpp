/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * IChol.cpp
 *
 *  Created on: Jan 29, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/IChol.hpp>

#include <math.h>
#include <iostream>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::SparseDataMatrix;

int main(int argc, char** argv) {
  //  const auto sweeps = [argc, argv]() {
  //    if (argc != 2) {
  //      printf("invalid amount of sweeps specified");
  //      exit(1);
  //    } else {
  //      try {
  //        return static_cast<size_t>(std::stoul(argv[1]));
  //      } catch (std::invalid_argument& e) {
  //        printf("amount of sweeps could not be parsed.\n%s\n", e.what());
  //        exit(1);
  //      }
  //    }
  //  }();
  //

  auto size = 5u;

  // initialize
  const std::vector<double> aSparse{2, 2, 6, 1, 6, 1, 6, 1, 4, 2, 16};
  const std::vector<size_t> colIdx{0, 0, 1, 0, 2, 2, 3, 0, 2, 3, 4};
  const std::vector<size_t> rowPtr{0, 1, 3, 5, 7};
  const std::vector<double> results{1.4142, 1.4142, 2.0000, 0.7071, 2.3452, 0.4264,
                                    2.4121, 0.7071, 1.4924, 0.5653, 3.5990};

  SparseDataMatrix matrix{size, size, aSparse, colIdx, rowPtr};

  DataMatrix A;
  SparseDataMatrix::toDataMatrix(matrix, A);
  std::cout << "before:\n" << A.toString() << "\n";

  // decomp:
  DataVector aNorm(size);
  sgpp::datadriven::IChol::normToUnitDiagonal(matrix, aNorm);
  SparseDataMatrix::toDataMatrix(matrix, A);
  std::cout << "normed:\n" << A.toString() << "\n";
  sgpp::datadriven::IChol::decompose(A, matrix, 1);
  // IChol::reaplyDiagonal(A, aNorm);

  SparseDataMatrix::toDataMatrix(matrix, A);
  std::cout << "after:\n" << A.toString() << "\n";

  //  const auto matSize = matrix.getNrows();
  //  // get the data vector
  //  auto& matData = matrix.getDataVector();
  //  // get the rows
  //  const auto& rowPtrs = matrix.getRowPtrVector();
  //  // get the cols
  //  const auto& colIndices = matrix.getColIndexVector();
  //
  //  const auto col = 3;
  //  const auto row = 4;
  //  const auto dataIter = 5;
  //
  //  auto s = matData[dataIter];
  //  std::cout << "s:" << s << "\n";
  //
  //  std::cout << rowPtrs[col] << " " << rowPtr[col + 1] - 1 << " " << rowPtrs[row] << " " <<
  //  dataIter
  //            << std::endl;
  //  auto upperFirst = colIndices.begin() + rowPtrs[col];
  //  const auto upperLast = colIndices.begin() + (rowPtr[col + 1] - 1);
  //  auto lowerFist = colIndices.begin() + rowPtrs[row];
  //  const auto lowerLast = colIndices.begin() + dataIter;
  //
  //  printf("upperFirst %d\nupperLast %d\nlowerFirst %d\nlowerLast %d\n", *upperFirst, *upperLast,
  //         *lowerFist, *lowerLast);
  //
  //  // sparse dot product by merging in O(n+m)
  //  while (lowerFist != lowerLast) {
  //    // if we're out of nonzeors in the upper row, then we're also done.
  //    if (upperFirst == upperLast) {
  //      std::cout << "done\n";
  //      break;
  //    }
  //    if (*upperFirst < *lowerFist) {
  //      printf("Mismatch: upper %d < %d lower: incrementing to ", *upperFirst, *lowerFist);
  //      ++upperFirst;
  //      std::cout << *upperFirst << std::endl;
  //    } else if (*lowerFist < *upperFirst) {
  //      printf("Mismatch: upper %d < %d lower: incrementing to ", *upperFirst, *lowerFist);
  //      ++lowerFist;
  //      std::cout << *lowerFist << std::endl;
  //    } else {
  //      printf("Match: upper %d == %d lower\n", *upperFirst, *lowerFist);
  //      printf("%f,%f\n", matData[upperFirst - colIndices.begin()],
  //             matData[lowerFist - colIndices.begin()]);
  //      s -= matData[upperFirst - colIndices.begin()] * matData[lowerFist - colIndices.begin()];
  //      std::cout << "s:" << s << "\n";
  //      ++upperFirst;
  //      ++lowerFist;
  //    }
  //  }
  //
  //  std::cout << "s:" << s << "\n";
  //
  //  if (row != col) {
  //    const auto index = row + 1 > matSize ? rowPtrs[row + 1] - 1 : matData.size() - 1;
  //    matData[dataIter] = s / matData[index];
  //  } else {
  //    matData[dataIter] = sqrt(s);
  //  }
  //
  //  SparseDataMatrix::toDataMatrix(matrix, A);
  //  std::cout << A.toString() << "\n";

  //  // test
  //  const auto aData = A.getDataVector();
  //  const auto aColIdx = A.getColIndexVector();
  //  const auto aPtr = A.getRowPtrVector();
  //
  //  DataVector aNorm{size};
  //  // norm:
  //  IChol::normToUnitDiagonal(A, aNorm);
  //
  //  DataMatrix m;
  //  SparseDataMatrix::toDataMatrix(A, m);
  //
  //  std::cout << m.toString() << "\n" << aNorm.toString();
}
