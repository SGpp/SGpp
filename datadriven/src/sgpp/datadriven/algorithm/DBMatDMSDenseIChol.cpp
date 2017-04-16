/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatDMSDenseIChol.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp>

namespace sgpp {
namespace datadriven {
using sgpp::base::DataMatrix;

// void sgpp::datadriven::DBMatDMSDenseIChol::choleskyUpdateLambda(
//    sgpp::base::DataMatrix& decompMatrix, double lambda_up) const {
//#pragma omp parallel
//#pragma omp single
//
//  DataMatrix systemMatrix{decompMatrix.getNrows(), decompMatrix.getNcols()};
//
//  { /*omp parallel region begin */
//// we need the full matrix again
//#pragma omp task
//    { /* omp task begin */
//
//    } /* omp task end */
//
//// meanwhile we can update the system matrix
//#pragma omp task
//    { /* omp task begin */
//#pragma omp simd
//      for (size_t i = 0; i < decompMatrix.getNrows(); i++) {
//        const auto val = decompMatrix.get(i, i);
//        decompMatrix.set(i, i, val + lambda_up);
//      }
//    } /* omp task end */
//  }   /*omp parallel region end */
//
//  // and run the decomposition again. Since out guess should be closer, we should converge much
//  // faster.
//  IChol::decompose(decompMatrix, decompMatrix, 2);
//}

} /* namespace datadriven */
} /* namespace sgpp */
