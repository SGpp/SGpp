// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <iomanip>
#include <string>

namespace sgpp {
namespace datadriven {

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const DBMatDensityConfiguration& config)
    : DBMatOffline(config) {
  this->lambda = config.lambda_;

  this->q_ortho_matrix_ = sgpp::base::DataMatrix(1, 1);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(1, 1);

  // dim_a = 0, indirectly tells the online object if build() or decompose() were performed
  this->dim_a = 0;
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
  // allocating subdiagonal and diagonal vectors of T
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

  // decomposed matrix: (lhs+lambda*I) = Q * T^{-1} * Q^t
  this->isDecomposed = true;
}

void DBMatOfflineOrthoAdapt::hessenberg_decomposition(gsl_vector* diag, gsl_vector* subdiag) {
  gsl_vector* tau = gsl_vector_alloc(dim_a - 1);
  gsl_matrix_view gsl_lhs = gsl_matrix_view_array(this->lhsMatrix.getPointer(), dim_a, dim_a);
  gsl_matrix_view gsl_q = gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), dim_a, dim_a);

  gsl_linalg_symmtd_decomp(&gsl_lhs.matrix, tau);
  gsl_linalg_symmtd_unpack(&gsl_lhs.matrix, tau, &gsl_q.matrix, diag, subdiag);
}

void DBMatOfflineOrthoAdapt::invert_symmetric_tridiag(gsl_vector* diag, gsl_vector* subdiag) {
  gsl_vector* e = gsl_vector_calloc(diag->size);  // calloc set all values to zero
  gsl_vector* x = gsl_vector_alloc(diag->size);   // target of solving

  for (size_t k = 0; k < this->t_tridiag_inv_matrix_.getNcols(); k++) {
    e->data[k] = 1;
    gsl_linalg_solve_symm_tridiag(diag, subdiag, e, x);
    for (size_t i = 0; i < this->t_tridiag_inv_matrix_.getNrows(); i++) {
      this->t_tridiag_inv_matrix_.set(k, i, x->data[i]);
    }
    e->data[k] = 0;
  }
}

void DBMatOfflineOrthoAdapt::store(const std::string& fileName) {
  DBMatOffline::store(fileName);

  FILE* outCFile = fopen(fileName.c_str(), "ab");
  if (!outCFile) {
    throw sgpp::base::algorithm_exception{"cannot open file for writing"};
  }

  // store q_ortho_matrix_
  gsl_matrix_view q_view =
      gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_fwrite(outCFile, &q_view.matrix);

  // store t_inv_tridiag_
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(this->t_tridiag_inv_matrix_.getPointer(), this->dim_a, this->dim_a);
  gsl_matrix_fwrite(outCFile, &t_inv_view.matrix);
  fclose(outCFile);
}
}  // namespace datadriven
}  // namespace sgpp
#endif /* USE_GSL */
