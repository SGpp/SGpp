// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>

#include <iomanip>
#include <iostream>
namespace sgpp {
namespace datadriven {

void DBMatDMSOrthoAdapt::solve(sgpp::base::DataMatrix& T_inv, sgpp::base::DataMatrix& Q,
                               sgpp::base::DataMatrix& B, sgpp::base::DataVector& b,
                               sgpp::base::DataVector& alpha) {
  // assert dimensions
  bool prior_refined = (B.getNcols() > 1);  // if B.getNcols <= 1, then not refined yet

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
  gsl_matrix_view b_vector_view_cut = gsl_matrix_view_array(b.getPointer(), Q.getNrows(), 1);
  gsl_matrix_view b_vector_view = gsl_matrix_view_array(b.getPointer(), b.getSize(), 1);
  gsl_matrix_view alpha_view_cut = gsl_matrix_view_array(alpha.getPointer(), Q.getNrows(), 1);
  gsl_matrix_view alpha_view = gsl_matrix_view_array(alpha.getPointer(), alpha.getSize(), 1);

  gsl_matrix* interim2 = gsl_matrix_alloc(Q.getNrows(), 1);

  // calculating: Q * T^{-1} * Q^t * b
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &q_view.matrix, &b_vector_view_cut.matrix, 0.0,
                 &alpha_view_cut.matrix);

  // calculating T^{-1} * Q^t * b
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &t_inv_view.matrix, &alpha_view_cut.matrix, 0.0,
                 interim2);

  // calculating Q * T^{-1} * Q^t * b
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &q_view.matrix, interim2, 0.0,
                 &alpha_view_cut.matrix);

  // if B should not be considered
  if (!prior_refined || B.getNcols() == Q.getNcols()) {
    if (interim2->size1 != alpha.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatDMSOrthoAdapt::solve: vector alpha does not match Q * T^{-1} * Q^t * b");
    }
  } else {
    // add the B*b term: alpha = alpha + B*b
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &b_matrix_view.matrix, &b_vector_view.matrix,
                   1.0, &alpha_view.matrix);  // gsl_blas allows to add to the target after mult.
  }

  // DEBUG: print alpha after solving
  // std::cout << "alpha after ortho_adapt_solve: \n";
  // for (size_t i = 0; i < alpha.getSize(); i++) {
  //   std::cout << alpha.get(i) << "   ";
  // }
  // std::cout << std::endl;
  gsl_matrix_free(interim2);
}
}  // namespace datadriven
}  // namespace sgpp
#endif /* USE_GSL */
