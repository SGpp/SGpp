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

#include <algorithm>
#include <list>
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
      DataMatrix matrix(lhsMatrix);
      IChol::decompose(matrix, lhsMatrix, 4);
    }
    isDecomposed = true;
  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
}

void DBMatOfflineDenseIChol::choleskyModification(size_t newPoints, std::list<size_t> deletedPoints,
                                                  double lambda) {
  std::cout << "dense ichol modification\n";
  if (newPoints > 0) {
    size_t gridSize = grid->getStorage().getSize();
    size_t gridDim = grid->getStorage().getDimension();

    // DataMatrix to collect vectors to append
    DataMatrix mat_refine(newPoints, gridSize);

    printf("mat size will be %d, %d\n", newPoints, gridSize);

    DataMatrix level(gridSize, gridDim);
    DataMatrix index(gridSize, gridDim);

    grid->getStorage().getLevelIndexArraysForEval(level, index);
    double lambda_conf = lambda;
// Loop to calculate all L2-products of added points based on the
// hat-function as basis function
#pragma omp parallel for schedule(guided)
    for (size_t j = gridSize - newPoints; j < gridSize; j++) {
      for (size_t i = 0; i <= j; i++) {
        double res = 1;
        for (size_t k = 0; k < gridDim; k++) {
          double lik = level.get(i, k);
          double ljk = level.get(j, k);
          double iik = index.get(i, k);
          double ijk = index.get(j, k);

          if (lik == ljk) {
            if (iik == ijk) {
              // Use formula for identical ansatz functions:
              res *= 2 / lik / 3;
            } else {
              // Different index, but same level => ansatz functions do not
              // overlap:
              res = 0.;
              break;
            }
          } else {
            if (std::max((iik - 1) / lik, (ijk - 1) / ljk) >=
                std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
              // Ansatz functions do not not overlap:
              res = 0.;
              break;
            } else {
              // Use formula for different overlapping ansatz functions:
              if (lik > ljk) {  // Phi_i_k is the "smaller" ansatz function
                double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
                double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
                temp_res *= ljk;
                temp_res = (1 - temp_res) / lik;
                res *= temp_res;
              } else {  // Phi_j_k is the "smaller" ansatz function
                double diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
                double temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
                temp_res *= lik;
                temp_res = (1 - temp_res) / ljk;
                res *= temp_res;
              }
            }
          }
        }
        // The new Rows/Cols are stored in mat_refine

        // add current lambda to lower diagonal elements of mat_refine
        if (i == j) {
          mat_refine.set(j - gridSize + newPoints, i, res + lambda_conf);
        } else {
          mat_refine.set(j - gridSize + newPoints, i, res);
        }
      }
    }

    std::cout << "mat_refine:\n" << mat_refine.toString() << "\n\n";
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
