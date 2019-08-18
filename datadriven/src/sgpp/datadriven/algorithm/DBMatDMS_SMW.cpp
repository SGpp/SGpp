// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDMS_SMW.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif /* USE_GSL */

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

void DBMatDMS_SMW::solve(sgpp::base::DataMatrix& A_inv, sgpp::base::DataMatrix& B,
                               sgpp::base::DataVector& b, sgpp::base::DataVector& alpha) {
#ifdef USE_GSL
  // assert dimensions
  bool prior_refined = (B.getNcols() > 1);  // if B.getNcols <= 1, then no refining yet

  if (prior_refined) {
    if (B.getNcols() != b.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SWM::solve: Matrix B and Vector b don't match for mult.");
    }
    if (B.getNcols() != alpha.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SWM::solve: Matrix B and vector alpha don't match for mult.");
    }
  }
  gsl_matrix_view A_inv_view =
      gsl_matrix_view_array(A_inv.getPointer(), A_inv.getNrows(), A_inv.getNcols());

  // here, b and alpha have to be adjusted to size of A_inv, to ignore 0 rows and columns
  gsl_vector_view b_cut_view = gsl_vector_view_array(b.getPointer(), A_inv.getNcols());
  gsl_vector_view alpha_cut_view = gsl_vector_view_array(alpha.getPointer(), A_inv.getNcols());

  // alpha = A^-1 * b
  gsl_blas_dgemv(CblasNoTrans, 1.0, &A_inv_view.matrix, &b_cut_view.vector, 0.0,
                 &alpha_cut_view.vector);

  if (A_inv.getNcols() < B.getNcols()) {
    gsl_matrix_view B_view = gsl_matrix_view_array(B.getPointer(), B.getNrows(), B.getNcols());
    gsl_vector_view b_view = gsl_vector_view_array(b.getPointer(), B.getNcols());
    gsl_vector_view alpha_view = gsl_vector_view_array(alpha.getPointer(), B.getNcols());

    // alpha = alpha + B*b
    gsl_blas_dgemv(CblasNoTrans, 1.0, &B_view.matrix, &b_view.vector, 1.0, &alpha_view.vector);
  }
#else
  throw sgpp::base::algorithm_exception("USE_GSL not set");
#endif /* USE_GSL */
}

void DBMatDMS_SMW::solveParallel(DataMatrixDistributed& T_inv, DataMatrixDistributed& Q,
                                       DataMatrixDistributed& B, DataVectorDistributed& b,
                                       DataVectorDistributed& alpha) {
#ifdef USE_SCALAPACK
  // assert dimensions
  bool prior_refined = (B.getGlobalCols() > 1);  // if B.getNcols <= 1, then no refining yet

  if (prior_refined) {
    if (B.getGlobalCols() != b.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SMW::solve: Matrix B and Vector b don't match for mult.");
    }
    if (B.getGlobalCols() != alpha.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SMW::solve: Matrix B and vector alpha don't match for mult.");
    }
  }

  auto processGrid = alpha.getProcessGrid();

  DataVectorDistributed tmp(processGrid, Q.getGlobalRows(), alpha.getBlockSize());

  /**
   * calculation of alpha = Q * T_inv * Q^t * b + B * b
   * note: T_inv and Q have only zeroes on indices reaching out the quadratic
   * size of  the original non refined grid, so the operations Q * T_inv * Q_t * b will be capped to
   * this size.
   */

  // calculating Q^t * b (result saved in part of alpha)
  DataMatrixDistributed::mult(Q, b, alpha, true);

  // calculating T^{-1} * Q^t * b (using part of alpha from previous line)
  DataMatrixDistributed::mult(T_inv, alpha, tmp, false);

  // calculating Q * T^{-1} * Q^t * b (result saved in part of alpha)
  DataMatrixDistributed::mult(Q, tmp, alpha, false);

  // if B should not be considered
  if (!prior_refined || B.getGlobalCols() == Q.getGlobalCols()) {
    if (tmp.getGlobalRows() != alpha.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SMW::solveParallel: vector alpha does not match Q * T^{-1} * Q^t * b");
    }
  } else {
    // add the B*b term: alpha = B * b + alpha (using beta parameter from pdgemv)
    DataMatrixDistributed::mult(B, b, alpha, false, 1.0, 1.0);
  }

  // DEBUG: print alpha after solving
  // alpha.printVector();
#else
  throw sgpp::base::algorithm_exception("USE_SCALAPACK not set");
#endif /* USE_SCALAPACK */
}
}  // namespace datadriven
}  // namespace sgpp
