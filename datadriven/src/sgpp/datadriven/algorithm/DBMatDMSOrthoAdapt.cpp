// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>

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

void DBMatDMSOrthoAdapt::solve(sgpp::base::DataMatrix& T_inv, sgpp::base::DataMatrix& Q,
                               sgpp::base::DataMatrix& B, sgpp::base::DataVector& b,
                               sgpp::base::DataVector& alpha) {
#ifdef USE_GSL
  // assert dimensions
  bool prior_refined = (B.getNcols() > 1);  // if B.getNcols <= 1, then no refining yet

  if (prior_refined) {
    if (B.getNcols() != b.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: Matrix B and Vector b don't match for mult.");
    }
    if (B.getNcols() != alpha.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: Matrix B and vector alpha don't match for mult.");
    }
  }

  /**
   * note: T_inv and Q have only zeroes on indices reaching out the quadratic
   * size of offline.dimA(), which is the size of the original non refined grid,
   * so the operations Q * T_inv * Q_t * b will be capped to size dimA()
   */

  // creating gsl_matrix_views to be able to use BLAS operations
  gsl_matrix_view q_view = gsl_matrix_view_array(Q.getPointer(), Q.getNrows(), Q.getNcols());
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(T_inv.getPointer(), T_inv.getNrows(), T_inv.getNcols());
  gsl_matrix_view b_matrix_view = gsl_matrix_view_array(B.getPointer(), B.getNrows(), B.getNcols());

  gsl_vector_view b_vector_view_cut = gsl_vector_view_array(b.getPointer(), Q.getNrows());
  gsl_vector_view b_vector_view = gsl_vector_view_array(b.getPointer(), b.getSize());
  gsl_vector_view alpha_view_cut = gsl_vector_view_array(alpha.getPointer(), Q.getNrows());
  gsl_vector_view alpha_view = gsl_vector_view_array(alpha.getPointer(), alpha.getSize());

  gsl_vector* interim2 = gsl_vector_alloc(Q.getNrows());

  // calculating Q^t * b
  gsl_blas_dgemv(CblasTrans, 1.0, &q_view.matrix, &b_vector_view_cut.vector, 0.0,
                 &alpha_view_cut.vector);

  // calculating T^{-1} * Q^t * b
  gsl_blas_dgemv(CblasNoTrans, 1.0, &t_inv_view.matrix, &alpha_view_cut.vector, 0.0, interim2);

  // calculating Q * T^{-1} * Q^t * b
  gsl_blas_dgemv(CblasNoTrans, 1.0, &q_view.matrix, interim2, 0.0, &alpha_view_cut.vector);

  // if B should not be considered
  if (!prior_refined || B.getNcols() == Q.getNcols()) {
    if (interim2->size != alpha.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: vector alpha does not match Q * T^{-1} * Q^t * b");
    }
  } else {
    // add the B*b term: alpha = alpha + B*b
    gsl_blas_dgemv(CblasNoTrans, 1.0, &b_matrix_view.matrix, &b_vector_view.vector, 1.0,
                   &alpha_view.vector);  // gsl_blas allows to add to the target after mult.
  }

  // DEBUG: print alpha after solving
  // std::cout << "alpha after ortho_adapt_solve: \n";
  // for (size_t i = 0; i < alpha.getSize(); i++) {
  //   std::cout << alpha.get(i) << "   ";
  // }
  // std::cout << std::endl;

  gsl_vector_free(interim2);
#else
  throw sgpp::base::algorithm_exception("USE_GSL not set");
#endif /* USE_GSL */
}

void DBMatDMSOrthoAdapt::solveParallel(DataMatrixDistributed& T_inv, DataMatrixDistributed& Q,
                                       DataMatrixDistributed& B, DataVectorDistributed& b,
                                       DataVectorDistributed& alpha) {
#ifdef USE_SCALAPACK
  // assert dimensions
  bool prior_refined = (B.getGlobalCols() > 1);  // if B.getNcols <= 1, then no refining yet

  if (prior_refined) {
    if (B.getGlobalCols() != b.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: Matrix B and Vector b don't match for mult.");
    }
    if (B.getGlobalCols() != alpha.getGlobalRows()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: Matrix B and vector alpha don't match for mult.");
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
          "In DBMatDMSOrthoAdapt::solveParallel: vector alpha does not match Q * T^{-1} * Q^t * b");
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
