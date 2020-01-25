// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>

//#include <ctime>
#include <iostream>

namespace sgpp {
namespace datadriven {

DBMatDMSBackSub::DBMatDMSBackSub() {}

DBMatDMSBackSub::~DBMatDMSBackSub() {}

void DBMatDMSBackSub::solve(sgpp::base::DataMatrix& DecompMatrix,
                            sgpp::base::DataVector& alpha,
                            sgpp::base::DataVector& b) {
  size_t resultSize = alpha.getSize();

  clock_t end;
  clock_t begin;

  //double elapsed_secs;

  //begin = clock();

  // Forward Substitution:
  sgpp::base::DataVector y(resultSize);
  for (size_t i = 0; i < resultSize; i++) {
    y[i] = b[i];
    for (size_t j = 0; j < i; j++) {
      y[i] -= DecompMatrix.get(i, j) * y[j];
    }
    // There is no need to divide by the diagonal element, because all
    // of them are 1 in L
  }

  // Backward Substitution:
  for (int i = static_cast<int>(resultSize) - 1; i >= 0; i--) {
    alpha[i] = y[i];
    for (size_t j = i + 1; j < resultSize; j++) {
      alpha[i] -= DecompMatrix.get(i, j) * alpha[j];
    }
    alpha[i] /= DecompMatrix.get(i, i);
  }

  //end = clock();
  //elapsed_secs = static_cast<double>(end - begin) / CLOCKS_PER_SEC;
  //std::cout << "Solve LU: " << elapsed_secs;
}

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
