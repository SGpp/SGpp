/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt.cpp
 *
 *  Created on: 01.08.2017
 *  Author: Dmitrij Boschko
 */

// #ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const DBMatDensityConfiguration& config)
    : DBMatOffline(config) {
  // todo: alternativ: berechne size ohne buildMatrix()
  this->buildMatrix();
  this->dim_a = this->getGrid().getStorage().getSize();
  this->lambda = config.lambda_;

  q_ortho_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
  t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
  diag_ = sgpp::base::DataVector(dim_a);
  subdiag_ = sgpp::base::DataVector(dim_a - 1);
}

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const std::string& fileName)
    : DBMatOffline(fileName) {
  // todo: Q and T_inv from file
}

DBMatOffline* DBMatOfflineOrthoAdapt::clone() { return new DBMatOfflineOrthoAdapt{*this}; }

bool DBMatOfflineOrthoAdapt::isRefineable() { return true; }

void DBMatOfflineOrthoAdapt::buildMatrix() {
  // build lhs matrix R
  DBMatOffline::buildMatrix();

  // set dimension of matrices Q, T_inv, etc.
  // this->dim_a = this->grid->getStorage().getSize(); // done in constructor

  isConstructed = true;
  /*
   * note: the regularization parameter lambda is added in the function
   * "invert_tridiag", right before the inverting computations start. This is
   * possible because R+lambda*I = Q*(T + lambda*I)*Q_t
   */
}

void DBMatOfflineOrthoAdapt::decomposeMatrix() {}

void DBMatOfflineOrthoAdapt::hessenberg_decomposition() {
  // todo:
  // currently this function uses gsl and relies on gsl_matrices, which means
  // it has to copy between gsl matrices and sgpp matrices, which is inefficient
  // gsl matrix view ?
  if (dim_a == 1) {
    this->q_ortho_matrix_.set(0, 0, 1.0);
    return;
  }
  gsl_vector* gsl_diag = gsl_vector_alloc(dim_a);
  gsl_vector* gsl_subdiag = gsl_vector_alloc(dim_a - 1);
  gsl_vector* tau = gsl_vector_alloc(dim_a - 1);
  gsl_matrix* gsl_lhs = gsl_matrix_alloc(dim_a, dim_a);
  gsl_matrix* gsl_q = gsl_matrix_alloc(dim_a, dim_a);

  for (size_t i = 0; i < dim_a; i++) {
    for (size_t j = 0; j < dim_a; j++) {
      gsl_matrix_set(gsl_lhs, i, j, this->lhsMatrix.get(i, j));
    }
  }

  gsl_linalg_symmtd_decomp(gsl_lhs, tau);
  gsl_linalg_symmtd_unpack(gsl_lhs, tau, gsl_q, gsl_diag, gsl_subdiag);

  // write computed values into members of class
  for (size_t i = 0; i < dim_a; i++) {
    (this->getDiag()).set(i, gsl_vector_get(gsl_diag, i));
    if (i < dim_a - 1) {
      this->subdiag_.set(i, gsl_vector_get(gsl_subdiag, i));
    }
    for (size_t j = 0; j < dim_a; j++) {
      this->q_ortho_matrix_.set(i, j, gsl_matrix_get(gsl_q, i, j));
    }
  }
  return;
}

// inverts a symmetric tridiag matrix
void DBMatOfflineOrthoAdapt::invert_tridiag() {
  size_t n = this->dim_a;

  // todo: copy constructor doesn't work?, therefore manual copy
  sgpp::base::DataVector superdiag(n, 0.0);
  for (size_t i = 0; i < n; i++) {
    superdiag.set(i, this->subdiag_.get(i));
  }

  // LR-decomposition
  for (size_t i = 0; i < n - 1; i++) {
    // l_i = l_i / a_i
    this->subdiag_.set(i, this->subdiag_.get(i) / this->diag_.get(i));
    // a_(i+1) = a_(i+1) - l_i*u_i
    this->diag_.set(i + 1, this->diag_.get(i + 1) - this->subdiag_.get(i) * superdiag.get(i));
  }

  // L*R*T_inv = Id, solve for columns of T_inv
  for (size_t i = 0; i < n; i++) {  // i-th column of T_inv
    // forward subst. L*y = e_i
    sgpp::base::DataVector y(n, 0.0);
    y.set(0, i == 0 ? 1.0 : 0.0);  // y_0 = (e_i)_0
    for (size_t j = 1; j < n; j++) {
      y.set(j, (i == j ? 1.0 : 0.0) - y.get(j - 1) * this->subdiag_.get(j - 1));
    }

    // backward subst. R*x = y, where x = columnf of T_inv
    this->t_tridiag_inv_matrix_.set(n - 1, i, y.get(n - 1) / this->diag_.get(n - 1));
    for (size_t j = n; j > 0; j--) {  // index-shift +1, because size_t >= 0
      double value = y.get(j - 1) - this->t_tridiag_inv_matrix_.get(j, i) * superdiag.get(j - 1);
      value /= this->diag_.get(j - 1);
      this->t_tridiag_inv_matrix_.set(j - 1, i, value);
    }
  }
}
}  // namespace datadriven
}  // namespace sgpp
// #endif /* USE_GSL */
