/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt.cpp
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <iomanip>
#include <string>

#ifdef USE_GSL

namespace sgpp {
namespace datadriven {

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const DBMatDensityConfiguration& config)
    : DBMatOffline(config) {
  this->dim_a = 1;
  this->lambda = config.lambda_;

  this->q_ortho_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
}

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const std::string& fileName)
    : DBMatOffline(fileName) {
  // grid already initialized in super constructor
  this->dim_a = this->getGrid().getStorage().getSize();
  this->lambda = this->config.lambda_;

  this->lhsMatrix = sgpp::base::DataMatrix(dim_a, dim_a);
  this->q_ortho_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);

  gsl_matrix_view lhs_view =
      gsl_matrix_view_array(this->lhsMatrix.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_view q_view =
      gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(this->t_tridiag_inv_matrix_.getPointer(), this->dim_a, this->dim_a);

  FILE* file = fopen(fileName.c_str(), "rb");
  if (!file) {
    throw sgpp::base::algorithm_exception{"Failed to open File"};
  }

  // seek end of first line
  char c = 0;
  while (c != '\n') {
    c = static_cast<char>(fgetc(file));
  }

  gsl_matrix_fread(file, &lhs_view.matrix);
  gsl_matrix_fread(file, &q_view.matrix);
  gsl_matrix_fread(file, &t_inv_view.matrix);

  fclose(file);

  this->isConstructed = true;
  this->isDecomposed = true;
}

DBMatOffline* DBMatOfflineOrthoAdapt::clone() { return new DBMatOfflineOrthoAdapt{*this}; }

bool DBMatOfflineOrthoAdapt::isRefineable() { return true; }

void DBMatOfflineOrthoAdapt::buildMatrix() {
  DBMatOffline::buildMatrix();
  this->dim_a = this->getGrid().getStorage().getSize();

  this->q_ortho_matrix_.resizeQuadratic(dim_a);
  this->t_tridiag_inv_matrix_.resizeQuadratic(dim_a);
}

void DBMatOfflineOrthoAdapt::decomposeMatrix() {
  // allocating sub-, super- and diagonal vectors of T
  gsl_vector* gsl_diag = gsl_vector_alloc(dim_a);
  gsl_vector* gsl_subdiag = gsl_vector_alloc(this->dim_a - 1);

  // decomposing: lhs = Q * T * Q^t
  this->hessenberg_decomposition(gsl_diag, gsl_subdiag);

  // adding configuration parameter lambda to diag before inverting T
  for (size_t i = 0; i < this->dim_a; i++) {
    gsl_diag->data[i] += this->lambda;
  }

  // inverting T+lambda*I, by solving L*R*x_i = e_i, for every i-th column x_i of T_inv
  this->invert_symmetric_tridiag(gsl_diag, gsl_subdiag);

  // decomposed matrix: (lhs+lambda*I) = Q * T_inv * Q^t
  this->isDecomposed = true;

  gsl_vector_free(gsl_diag);
  gsl_vector_free(gsl_subdiag);
}

void DBMatOfflineOrthoAdapt::hessenberg_decomposition(gsl_vector* diag, gsl_vector* subdiag) {
  gsl_vector* tau = gsl_vector_alloc(dim_a - 1);
  gsl_matrix_view gsl_lhs = gsl_matrix_view_array(this->lhsMatrix.getPointer(), dim_a, dim_a);
  gsl_matrix_view gsl_q = gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), dim_a, dim_a);

  gsl_linalg_symmtd_decomp(&gsl_lhs.matrix, tau);
  gsl_linalg_symmtd_unpack(&gsl_lhs.matrix, tau, &gsl_q.matrix, diag, subdiag);

  gsl_vector_free(tau);
  return;
}

void DBMatOfflineOrthoAdapt::invert_symmetric_tridiag(gsl_vector* diag, gsl_vector* subdiag) {
  /*
   * note: instead gsl_vector_get/set, gsl_vector->data[index] is used
   * this eliminates unnecessary range checking
   */
  const size_t n = this->dim_a;       // size of quadratic matrix to invert
  double* superdiag = new double[n];  // superdiagonal of T

  // LR-decomposition:
  // constructing T = LR and superdiag out of the vectors diag and subdiag
  for (size_t i = 0; i < n - 1; i++) {
    // copy subdiag into superdiag for later uses of old values
    superdiag[i] = subdiag->data[i];
    // l_i = l_i / a_i
    subdiag->data[i] = superdiag[i] / diag->data[i];
    // a_{i+1} = a_{i+1} - l_i*u_i
    diag->data[i + 1] = diag->data[i + 1] - subdiag->data[i] * superdiag[i];
  }
  // L*R*T_inv = Id, solve for columns x_i of T_inv
  for (size_t i = 0; i < n; i++) {
    //
    // forward subst. L*y = e_i
    double* y = new double[n];
    // y_0 = {e_i}_0
    y[0] = (i == 0 ? 1.0 : 0.0);
    // y_j = {e_i}_j - l_{j-1} * y_{j-1}
    for (size_t j = 1; j < n; j++) {
      y[j] = (i == j ? 1.0 : 0.0) - subdiag->data[j - 1] * y[j - 1];
    }

    // backward subst. R*x_i = y
    // {x_i}_{n-1} = y_{n-1} / a_{n-1}
    this->t_tridiag_inv_matrix_.set(n - 1, i, y[n - 1] / diag->data[n - 1]);
    // {x_i}_j = (y_j - u_j * x_{j+1} / a_j
    for (size_t j = n - 2; j >= i; j--) {  // j >= i, because symmetric
      double value =
          (y[j] - superdiag[j] * this->t_tridiag_inv_matrix_.get(j + 1, i)) / diag->data[j];
      this->t_tridiag_inv_matrix_.set(j, i, value);
      this->t_tridiag_inv_matrix_.set(i, j, value);
      // j >= i >= 0, because size_t, but i++ makes this only relevant for i==0
      if (j == 0) break;
    }
    delete[] y;
  }
  delete[] superdiag;
  return;
}

void DBMatOfflineOrthoAdapt::store(const std::string& fileName) {
  DBMatOffline::store(fileName);

  // #ifdef USE_GSL
  FILE* outCFile = fopen(fileName.c_str(), "ab");
  if (!outCFile) {
    throw sgpp::base::algorithm_exception{"cannot open file for writing"};
  }

  // store q_ortho_matrix_
  gsl_matrix_view q_view =
      gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_fwrite(outCFile, &q_view.matrix);

  // store t_inv_tridiag
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(this->t_tridiag_inv_matrix_.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_fwrite(outCFile, &t_inv_view.matrix);
  fclose(outCFile);
  // #else
  // throw base::not_implemented_exception("built withot GSL");
  // #endif /* USE_GSL */
}
}  // namespace datadriven
}  // namespace sgpp
#else
throw sgpp::base::algorithm_exception("USE_GSL is not set to true");
#endif /* USE_GSL */
