// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDMS_SMW.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

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

  // if B should not be considered (i.e. prior to refinement)
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

void DBMatDMS_SMW::solveParallel(DataMatrixDistributed& A_inv, DataMatrixDistributed& B,
                                 DataVectorDistributed& b, DataVectorDistributed& alpha) {
#ifdef USE_SCALAPACK
  // assert dimensions
  bool prior_refined = (B.getGlobalCols() > 1);  // if B.getNcols <= 1, then no refining yet

  if (prior_refined) {
    if (B.getGlobalRows() != b.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SWM::solveParallel: global size of matrix B and Vector b don't match for "
          "mult.");
    }
    if (B.getGlobalCols() != alpha.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMS_SWM::solveParallel: global size of matrix B and vector alpha don't match "
          "for mult.");
    }
  }

  auto processGrid = alpha.getProcessGrid();

  DataVectorDistributed tmp(processGrid, A_inv.getGlobalRows(), alpha.getBlockSize());

  // alpha = A^-1 * b (result saved in part of alpha)
  DataMatrixDistributed::mult(A_inv, b, alpha, true);

  // if B should not be considered (i.e. prior to refinement)
  if (A_inv.getGlobalCols() < B.getGlobalCols()) {
    // add the B*b term: alpha = B * b + alpha (using beta parameter from pdgemv)
    DataMatrixDistributed::mult(B, b, alpha, false, 1.0, 1.0);
  }
#else
  throw sgpp::base::algorithm_exception("USE_SCALAPACK not set");
#endif /* USE_SCALAPACK */
}
}  // namespace datadriven
}  // namespace sgpp
