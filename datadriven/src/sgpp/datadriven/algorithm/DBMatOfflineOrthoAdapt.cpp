// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif /* USE_GSL */

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::base::AdpativityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig)
    : DBMatOffline(gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig) {
  this->lambda = regularizationConfig.lambda_;

  this->q_ortho_matrix_ = sgpp::base::DataMatrix(1, 1);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(1, 1);

  // dim_a = 0, indirectly tells the online object if build() or decompose() were performed
  this->dim_a = 0;
}

DBMatOfflineOrthoAdapt::DBMatOfflineOrthoAdapt(const std::string& fileName)
    : DBMatOffline(fileName) {
  // grid already initialized in super constructor
  this->dim_a = this->getGrid().getStorage().getSize();
  this->lambda = this->regularizationConfig.lambda_;

  this->lhsMatrix = sgpp::base::DataMatrix(dim_a, dim_a);
  this->q_ortho_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
  this->t_tridiag_inv_matrix_ = sgpp::base::DataMatrix(dim_a, dim_a);
#ifdef USE_GSL
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
#else
  throw sgpp::base::algorithm_exception("USE_GSL has to be set");
#endif /* USE_GSL */
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
#ifdef USE_GSL
  // allocating subdiagonal and diagonal vectors of T
  sgpp::base::DataVector diag(dim_a);
  sgpp::base::DataVector subdiag(this->dim_a - 1);

  // decomposing: lhs = Q * T * Q^t
  this->hessenberg_decomposition(diag, subdiag);

  // adding configuration parameter lambda to diag before inverting T
  for (size_t i = 0; i < this->dim_a; i++) {
    diag.set(i, diag.get(i) + this->lambda);
  }

  // inverting T+lambda*I, by solving L*R*x_i = e_i, for every i-th column x_i of T_inv
  this->invert_symmetric_tridiag(diag, subdiag);

  // decomposed matrix: (lhs+lambda*I) = Q * T^{-1} * Q^t
  this->isDecomposed = true;
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::hessenberg_decomposition(sgpp::base::DataVector& diag,
                                                      sgpp::base::DataVector& subdiag) {
#ifdef USE_GSL
  gsl_vector* tau = gsl_vector_alloc(dim_a - 1);
  gsl_matrix_view gsl_lhs = gsl_matrix_view_array(this->lhsMatrix.getPointer(), dim_a, dim_a);
  gsl_matrix_view gsl_q = gsl_matrix_view_array(this->q_ortho_matrix_.getPointer(), dim_a, dim_a);
  gsl_vector_view gsl_diag = gsl_vector_view_array(diag.getPointer(), dim_a);
  gsl_vector_view gsl_subdiag = gsl_vector_view_array(subdiag.getPointer(), dim_a - 1);

  // does the decomposition
  gsl_linalg_symmtd_decomp(&gsl_lhs.matrix, tau);

  // unpacks information out of matrix to explicitly create Q, and T
  gsl_linalg_symmtd_unpack(&gsl_lhs.matrix, tau, &gsl_q.matrix, &gsl_diag.vector,
                           &gsl_subdiag.vector);

  gsl_vector_free(tau);
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::invert_symmetric_tridiag(sgpp::base::DataVector& diag,
                                                      sgpp::base::DataVector& subdiag) {
#ifdef USE_GSL
  gsl_vector* e = gsl_vector_calloc(diag.getSize());  // calloc sets all values to zero
  gsl_vector* x = gsl_vector_alloc(diag.getSize());   // target of solving

  gsl_vector_view gsl_diag = gsl_vector_view_array(diag.getPointer(), dim_a);
  gsl_vector_view gsl_subdiag = gsl_vector_view_array(subdiag.getPointer(), dim_a - 1);

  // loops columns of T_inv
  for (size_t k = 0; k < this->t_tridiag_inv_matrix_.getNcols(); k++) {
    e->data[k] = 1;
    gsl_linalg_solve_symm_tridiag(&gsl_diag.vector, &gsl_subdiag.vector, e, x);
    for (size_t i = 0; i < this->t_tridiag_inv_matrix_.getNrows(); i++) {
      this->t_tridiag_inv_matrix_.set(k, i, x->data[i]);
    }
    e->data[k] = 0;
  }

  gsl_vector_free(e);
  gsl_vector_free(x);
#endif /* USE_GSL */
}

void DBMatOfflineOrthoAdapt::store(const std::string& fileName) {
#ifdef USE_GSL
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
#endif /* USE_GSL */
}
}  // namespace datadriven
}  // namespace sgpp
