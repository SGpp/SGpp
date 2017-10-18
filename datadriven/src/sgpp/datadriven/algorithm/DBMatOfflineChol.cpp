/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineChol.cpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>
#endif /* USE_GSL */

#include <algorithm>
#include <chrono>
#include <list>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::algorithm_exception;

DBMatOfflineChol::DBMatOfflineChol(const DBMatDensityConfiguration& oc) : DBMatOfflineGE(oc) {}

DBMatOfflineChol::DBMatOfflineChol(const std::string& fileName) : DBMatOfflineGE{fileName} {}

DBMatOffline* DBMatOfflineChol::clone() { return new DBMatOfflineChol{*this}; }

bool DBMatOfflineChol::isRefineable() { return true; }

void DBMatOfflineChol::decomposeMatrix() {
#ifdef USE_GSL
  if (isConstructed) {
    if (isDecomposed) {
      // Already decomposed => Do nothing
      return;
    } else {
      auto begin = std::chrono::high_resolution_clock::now();

      size_t n = lhsMatrix.getNrows();
      gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                                n);  // Create GSL matrix view for decomposition
      // Perform Cholesky decomposition
      gsl_linalg_cholesky_decomp(&m.matrix);

      // Isolate lower triangular matrix
      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          if (i < j) {
            lhsMatrix.set(i, j, 0);
          }
        }
      }
      isDecomposed = true;
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Chol decomp took "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << "ms" << std::endl;
    }
  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
#else
  throw algorithm_exception("built withot GSL");
#endif /*USE_GSL*/
}

void DBMatOfflineChol::choleskyModification(size_t newPoints, std::list<size_t> deletedPoints,
                                            double lambda) {
#ifdef USE_GSL

  // Start coarsening
  // If list 'deletedPoints' is not empty, grid points got removed
  if (deletedPoints.size() > 0) {
    size_t new_size = grid->getSize();
    size_t old_size = new_size - newPoints + deletedPoints.size();

    // 'c' is the threshold to decide what kind of Cholesky modification should
    // be performed to take
    // account for the removed grid points
    double c = static_cast<double>(old_size) / 10;
    // Indicates which permutation to perform
    bool start1 = false;
    size_t index_coarse = 0;

    // If index of removed grid point <= c-> Permutation towards first
    // row/column -> job = 1
    size_t coarseCount_1 = 0;
    // If index of removed grid point > c -> Permutation towards last row/column
    // -> job = 2
    size_t coarseCount_2 = 0;

    for (std::list<size_t>::reverse_iterator it = deletedPoints.rbegin();
         it != deletedPoints.rend(); it++) {
      // Indicates current row/column to permute
      index_coarse = *it + 1 + coarseCount_1;

      if (static_cast<double>(index_coarse) > c && start1 == false) {
        // Is accessed if index is larger than 'c'
        choleskyPermutation(index_coarse, old_size - coarseCount_2, 2);
        coarseCount_2++;
      } else if (start1 == false) {
        // Is accessed the first 'index_coarse' is less than 'c'
        start1 = true;
        // Delete last 'coarseCount_2' rows/columns of current Cholesky factor
        lhsMatrix.resizeQuadratic(old_size - coarseCount_2);
      }

      if (start1 == true) {
        // Is accessed if 'index_coarse' is less than 'c'
        choleskyPermutation(1 + coarseCount_1, index_coarse, 1);
        coarseCount_1++;
      }
    }

    if (coarseCount_1 > 0) {
      // If some indices have been less than 'c'
      DataMatrix update_matrix(lhsMatrix.getNrows(), lhsMatrix.getNcols());
      update_matrix.copyFrom(lhsMatrix);
      // Resize copy of current Cholesky factor to
      // receive a matrix of rank one update vectors
      update_matrix.resizeToSubMatrix(coarseCount_1 + 1, 1, lhsMatrix.getNrows(), coarseCount_1);
      // Resize current Cholesky factor to required submatrix
      // for necessary rank one updates
      lhsMatrix.resizeToSubMatrix(coarseCount_1 + 1, coarseCount_1 + 1, lhsMatrix.getNrows(),
                                  lhsMatrix.getNrows());
      DataVector temp_col(update_matrix.getNrows());

      // 'coarseCount_1' many rank one updates based on the columns of
      // 'update_matrix' are performed
      DBMatDMSChol cholsolver;
      for (size_t i = 0; i < coarseCount_1; i++) {
        update_matrix.getColumn(i, temp_col);
        cholsolver.choleskyUpdate(lhsMatrix, temp_col, false);
      }
    } else {
      // If no indices have been less than 'c'
      lhsMatrix.resizeQuadratic(old_size - coarseCount_2);
    }
  }

  // Start refinement
  if (newPoints > 0) {
    size_t gridSize = grid->getStorage().getSize();
    size_t gridDim = grid->getStorage().getDimension();

    // DataMatrix to collect vectors to append
    DataMatrix mat_refine(gridSize, newPoints);

    DataMatrix level(gridSize, gridDim);
    DataMatrix index(gridSize, gridDim);

    grid->getStorage().getLevelIndexArraysForEval(level, index);
    double lambda_conf = lambda;
    // Loop to calculate all L2-products of added points based on the
    // hat-function as basis function
    for (size_t i = 0; i < gridSize; i++) {
      for (size_t j = gridSize - newPoints; j < gridSize; j++) {
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
          mat_refine.set(i, j - gridSize + newPoints, res + lambda_conf);
        } else {
          mat_refine.set(i, j - gridSize + newPoints, res);
        }
      }
    }

    // std::cout << "mat_refine:\n" << mat_refine.toString() << "\n\n";

    // Resize Cholesky factor to new 'gridSize' before 'choleskyAddPoint' is
    // applied
    // in order to save runtime
    this->lhsMatrix.resizeQuadratic(gridSize);

    // Now its time to call 'choleskyAddPoint''countNewGridPoints' often
    DataVector temp_col = DataVector(gridSize);
    for (size_t j = gridSize - newPoints; j < gridSize; j++) {
      temp_col.resizeZero(gridSize);
      mat_refine.getColumn(j - gridSize + newPoints, temp_col);
      temp_col.resizeZero(j + 1);
      choleskyAddPoint(temp_col, j);
    }
  }
#else
  throw algorithm_exception("built withot GSL");
#endif /*USE_GSL*/
}


void DBMatOfflineChol::choleskyAddPoint(DataVector& newCol, size_t size) {
#ifdef USE_GSL
  if (!isDecomposed) {
    throw algorithm_exception("Matrix was not decomposed, yet!");
  }

  DataMatrix& mat = lhsMatrix;
  // Size of provided memory for Cholesky factor,
  // because the allocations take place in 'choleskyModifications'
  size_t size_full = mat.getNrows();
  // Size of Cholesky factor after adding 'newCol'
  size_t size_up = newCol.getSize();

  if (size_up != (size + 1)) {
    throw algorithm_exception(
        "Size of update vector newCol needs to be 1 dim larger than the "
        "underlying decomposed matrix!");
  }

  // Create GSL matrix view for update procedures
  gsl_matrix_view m_full = gsl_matrix_view_array(mat.getPointer(), size_full, size_full);
  // Access submatrx since mat_ has already expanded to save alocations
  // procedures
  gsl_matrix_view m = gsl_matrix_submatrix(&m_full.matrix, 0, 0, size, size);

  gsl_vector_view vvec = gsl_vector_view_array(newCol.getPointer(), size_up);
  gsl_vector* wkvec_a = gsl_vector_calloc(size);

  // Extract newCol(0:size-1) = a
  gsl_vector_view c = gsl_vector_subvector(&vvec.vector, 0, size);

  // Solve system a = Lc
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m.matrix, &c.vector);
  gsl_blas_dcopy(&c.vector, wkvec_a);

  double phi;
  gsl_blas_ddot(wkvec_a, &c.vector, &phi);

  // Compute d = sqrt(newCol[last] - phi)
  double last = gsl_vector_get(&vvec.vector, size_up - 1);
  phi = last - phi;
  // Ensure 'phi' is larger 0
  if (phi <= 0) {
    throw algorithm_exception("Resulting matrix is at least not numerical positive definite!");
  } else {
    phi = sqrt(phi);
  }

  // Modify Choleskyfactor m (L)

  // DataVector full of zeros
  DataVector zeros(size_full, 0.0);

  // Add 'newCol' and 'zeros' to current Cholesky factor 'lhsMatrix_'
  mat.setColumn(size, zeros);
  newCol.set(size_up - 1, phi);
  newCol.resizeZero(size_full);
  mat.setRow(size, newCol);

  gsl_vector_free(wkvec_a);
#else
  throw algorithm_exception("built withot GSL");
#endif /*USE_GSL*/
}

void DBMatOfflineChol::choleskyPermutation(size_t k, size_t l, size_t job) {
#ifdef USE_GSL
  if (!isDecomposed) {
    throw algorithm_exception("Matrix was not decomposed, yet!");
  }

  DataMatrix& mat = lhsMatrix;
  size_t size = mat.getNrows();

  if (k > l) {
    throw algorithm_exception("l needs to be larger than k");
  } else if (l > size) {
    throw algorithm_exception("l needs to be smaller than the column size of the matrix");
  } else if (l == k) {
    return;
  }

  // Determine amount of necessary Givens rotations and rows to swap
  size_t count_perm = l - k;
  // Create GSL matrix view for update procedures
  gsl_matrix_view m = gsl_matrix_view_array(mat.getPointer(), size, size);

  // Define and declare Workingvector, Cosine- and Sinevector
  gsl_vector* svec = gsl_vector_calloc(count_perm);
  gsl_vector* cvec = gsl_vector_calloc(count_perm);

  double* givens_zero;
  double* tbuff = m.matrix.data;

  if (job == 2) {
    // Permute upper triangular - job = 2        => left circular shift
    // 1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ...,
    // l-1,l,k,l+1,..., size
    for (size_t i = k - 1; i < l - 1; i++) {
      gsl_matrix_swap_rows(&m.matrix, i, i + 1);
    }

    tbuff += (k - 1) * size + (k - 1);
    for (size_t j = 1; j <= count_perm; j++) {
      givens_zero = gsl_matrix_ptr(&m.matrix, k + j - 2, k + j - 1);
      // Givensrotation in (k+i-1,k+i)-Plane
      gsl_blas_drotg(tbuff, givens_zero, cvec->data + j - 1, svec->data + j - 1);
      // Access columns to modify via Givens rotation
      gsl_vector_view diag_sub =
          gsl_matrix_subcolumn(&m.matrix, k + j - 2, k + j - 1, size - k - j + 1);
      gsl_vector_view low_diag_sub =
          gsl_matrix_subcolumn(&m.matrix, k + j - 1, k + j - 1, size - k - j + 1);
      gsl_matrix_set(&m.matrix, k + j - 2, k + j - 1, 0.0);
      // Apply Givens rotation
      gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j - 1], svec->data[j - 1]);
      tbuff += (size + 1);
    }
  } else if (job == 1) {
    // Permute upper triangular - job = 1       => right circular shift
    // 1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ...,
    // l-1,l+1,...size
    for (size_t i = l - 1; i > k - 1; i--) {
      gsl_matrix_swap_rows(&m.matrix, i, i - 1);
    }

    tbuff += (k - 1) * size + (l - 2);
    for (size_t j = 1; j <= count_perm; j++) {
      givens_zero = gsl_matrix_ptr(&m.matrix, k - 1, l - j);
      // Givensrotation in (l-i,l-i+1)-Plane
      gsl_blas_drotg(tbuff, givens_zero, cvec->data + j - 1, svec->data + j - 1);
      // Access columns to modify via Givens rotation
      gsl_vector_view diag_sub = gsl_matrix_subcolumn(&m.matrix, l - j - 1, l - j, size - l + j);
      gsl_vector_view low_diag_sub = gsl_matrix_subcolumn(&m.matrix, l - j, l - j, size - l + j);
      gsl_matrix_set(&m.matrix, k - 1, l - j, 0.0);
      // Apply Givens rotation
      gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j - 1], svec->data[j - 1]);
      tbuff -= 1;
    }
  }

  gsl_vector_free(svec);
  gsl_vector_free(cvec);
#else
  throw algorithm_exception("built withot GSL");
#endif /*USE_GSL*/
}
} /* namespace datadriven */
} /* namespace sgpp */
