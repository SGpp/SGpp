// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>

#include <algorithm>
#include <chrono>
#include <list>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineDenseIChol::DBMatOfflineDenseIChol() : DBMatOfflineChol() {}

DBMatOfflineDenseIChol::DBMatOfflineDenseIChol(const std::string& fileName)
    : DBMatOfflineChol{fileName} {}

DBMatOffline* DBMatOfflineDenseIChol::clone() const { return new DBMatOfflineDenseIChol{*this}; }

void DBMatOfflineDenseIChol::decomposeMatrix(
    const RegularizationConfiguration& regularizationConfig,
    const DensityEstimationConfiguration& densityEstimationConfig) {
  if (isConstructed) {
    if (isDecomposed) {
      return;
    } else {
      //  auto begin = std::chrono::high_resolution_clock::now();

      DataMatrix tmpMatrix(lhsMatrix.getNrows(), lhsMatrix.getNcols());

// only copy lower triangular matrix
#pragma omp parallel for schedule(guided)
      for (size_t i = 0u; i < tmpMatrix.getNrows(); i++) {
#pragma omp simd
        for (size_t j = 0u; j <= i; j++) {
          tmpMatrix.set(i, j, lhsMatrix.get(i, j));
        }
      }

      ichol(tmpMatrix, lhsMatrix, densityEstimationConfig.iCholSweepsDecompose_);
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

void DBMatOfflineDenseIChol::choleskyModification(
    Grid& grid, datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    size_t newPoints, std::list<size_t> deletedPoints, double lambda, double c) {
  if (newPoints > 0) {
    //    auto begin = std::chrono::high_resolution_clock::now();

    size_t gridSize = grid.getSize();
    size_t gridDim = grid.getDimension();

    DataMatrix tmpMatrix(std::move(lhsMatrix));

    lhsMatrix = DataMatrix(gridSize, gridSize);

// only copy lower triangular matrix
#pragma omp parallel for schedule(guided)
    for (size_t i = 0u; i < tmpMatrix.getNrows(); i++) {
#pragma omp simd
      for (size_t j = 0u; j <= i; j++) {
        lhsMatrix.set(i, j, tmpMatrix.get(i, j));
      }
    }

    // lhsMatrix.resizeQuadratic(gridSize);

    // DataMatrix to collect newly added points to avoid a full copy of the matrix.
    DataMatrix matRefine(newPoints, gridSize);

    // printf("mat size will be %d, %d\n", newPoints, gridSize);

    DataMatrix level(gridSize, gridDim);
    DataMatrix index(gridSize, gridDim);

    grid.getStorage().getLevelIndexArraysForEval(level, index);
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
          lhsMatrix.set(j, i, res + lambda_conf);
          matRefine.set(j - (gridSize - newPoints), i, res + lambda_conf);
        } else {
          lhsMatrix.set(j, i, res);
          matRefine.set(j - (gridSize - newPoints), i, res);
        }
      }
    }

    ichol(matRefine, lhsMatrix, densityEstimationConfig.iCholSweepsRefine_, (gridSize - newPoints));

    //    auto end = std::chrono::high_resolution_clock::now();
    //    std::cout << "IChol refinement took "
    //              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<
    //              "ms"
    //              << std::endl;
  }
}

void sgpp::datadriven::DBMatOfflineDenseIChol::ichol(const DataMatrix& matrix, DataMatrix& result,
                                                     size_t sweeps, size_t startRow) {
  if (startRow > matrix.getSize()) {
    throw algorithm_exception{"Start row is larger then the matrix size"};
  }

// for all sweeps
#pragma omp parallel
  { /* omp parallel */
    for (size_t sweep = 0u; sweep < sweeps; sweep++) {
      // for each row
      for (size_t i = startRow; i < result.getNrows(); i++) {
// in each column until diagonal element
#pragma omp for schedule(guided) nowait
        for (size_t j = 0u; j <= i; j++) {
          // calculate sum;
          auto s = matrix.get(i - startRow, j);
          if (s > 0.0) {
#pragma omp simd
            for (size_t k = 0u; k < j; k++) {
              s -= result.get(i, k) * result.get(j, k);
            }
            if (i != j) {
              result.set(i, j, s / result.get(j, j));
            } else {
              result.set(i, i, sqrt(s));
            }
          }
        }
      }
    }
  }
} /* omp parallel */

sgpp::datadriven::MatrixDecompositionType DBMatOfflineDenseIChol::getDecompositionType() {
  return sgpp::datadriven::MatrixDecompositionType::DenseIchol;
}
} /* namespace datadriven */
} /* namespace sgpp */
