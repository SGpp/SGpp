// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>
#endif /* USE_GSL */

#include <math.h>
#include <ctime>
#include <iostream>

namespace sgpp {
namespace datadriven {

void DBMatDMSChol::solve(sgpp::base::DataMatrix& decompMatrix, sgpp::base::DataVector& alpha,
                         const sgpp::base::DataVector& b, double lambda_old,
                         double lambda_new) const {
  size_t size = decompMatrix.getNcols();
  // Performe Update based on Cholesky - afterwards perform n (GridPoints) many
  // rank-One-updates

  double lambda_up = lambda_new - lambda_old;

  // If regularization paramter is changed enter
  if (lambda_up != 0.0) {
    choleskyUpdateLambda(decompMatrix, lambda_up);
  }

  // Solve (R + lambda * I)alpha = b to obtain density declaring coefficents
  // alpha.

  // Forward Substitution:
  sgpp::base::DataVector y(size);
  choleskyForwardSolve(decompMatrix, b, y);

  // Backward Substitution:
  choleskyBackwardSolve(decompMatrix, y, alpha);

  std::cout << alpha.toString() << std::endl;
}

void DBMatDMSChol::solveParallel(DataMatrixDistributed& decompMatrix, DataVectorDistributed& x,
                                 double lambda_old, double lambda_new) const {
#ifdef USE_SCALAPACK

  // Performe Update based on Cholesky - afterwards perform n (GridPoints) many
  // rank-One-updates
  double lambda_up = lambda_new - lambda_old;

  // If regularization paramter is changed enter
  if (lambda_up != 0.0) {
    // current implementation: gather and update decomposition on master process
    DataMatrix decompMatrixLocal = decompMatrix.toLocalDataMatrix();

    auto processGrid = decompMatrix.getProcessGrid();
    if (processGrid->getCurrentRow() == 0 && processGrid->getCurrentColumn() == 0) {
      choleskyUpdateLambda(decompMatrixLocal, lambda_up);
    }

    decompMatrix.distribute(decompMatrixLocal.data(), 0, 0);
  }

  // Solve (R + lambda * I)alpha = b to obtain density declaring coefficents
  // alpha.

  DataMatrixDistributed::solveCholesky(decompMatrix, x);

  // DEBUG: print alpha after solving
  // x.printVector();
#else
  throw sgpp::base::algorithm_exception("build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

// Implement cholesky Update for given Decomposition and update vector
void DBMatDMSChol::choleskyUpdate(sgpp::base::DataMatrix& decompMatrix,
                                  const sgpp::base::DataVector& update, bool do_cv) const {
#ifdef USE_GSL
  // int i;

  size_t size = decompMatrix.getNrows();

  if (update.getSize() != size) {
    throw sgpp::base::data_exception(
        "choleskyUpdate::Size of DecomposedMatrix and updateVector don´t "
        "match...");
  }

  // Create GSL matrix view for update procedures
  gsl_matrix_view m = gsl_matrix_view_array(decompMatrix.getPointer(), size, size);
  gsl_vector_const_view vvec = gsl_vector_const_view_array(update.getPointer(), update.getSize());

  // Define and declare Workingvector, Cosine- and Sinevector
  gsl_vector* wkvec = gsl_vector_calloc(update.getSize());
  gsl_vector* svec = gsl_vector_calloc(update.getSize());
  gsl_vector* cvec = gsl_vector_calloc(update.getSize());

  // Generate Givens rotations, update L
  // Copy Values of updateVector into WorkingVector
  gsl_blas_dcopy(&vvec.vector, wkvec);

  double temp;
  double* tbuff = m.matrix.data;
  bool first_notZero = true;
  for (size_t i = 0; i < size - 1; i++) {
    if (*tbuff == 0.0 && wkvec->data[i] == 0.0) {
      throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
    } else if (wkvec->data[i] == 0.0 && do_cv == true) {
      // If cross validation is applied the first (n-1) entries in
      // the n-th step are zero and therefore can be skipped.
      tbuff += (size + 1);
      continue;
    } else if (wkvec->data[i] == 0.0 && first_notZero == true) {
      first_notZero = false;
      tbuff += (size + 1);
      continue;
    }
    // Determine givens rotation
    gsl_blas_drotg(tbuff, wkvec->data + i, cvec->data + i, svec->data + i);
    if ((temp = *tbuff) < 0.0) {
      *tbuff = -temp;
      cvec->data[i] = -cvec->data[i];
      svec->data[i] = -svec->data[i];
    } else if (temp == 0.0) {
      throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
    }

    // Access columns to modify via Givens rotations
    gsl_vector_view mat_sub = gsl_matrix_subcolumn(&m.matrix, i, i + 1, size - i - 1);
    gsl_vector_view wkvec_sub = gsl_vector_subvector(wkvec, i + 1, size - i - 1);

    // Allpy Givens rotation to mat_sub und wkvec_sub
    for (size_t j = 0; j < size - i - 1; j++) {
      double x = mat_sub.vector.data[mat_sub.vector.stride * j];
      double y = wkvec_sub.vector.data[j];
      mat_sub.vector.data[mat_sub.vector.stride * j] = cvec->data[i] * x + svec->data[i] * y;
      wkvec_sub.vector.data[j] = -svec->data[i] * x + cvec->data[i] * y;
    }
    tbuff += (size + 1);
  }

  size_t i_N = size - 1;
  // Apply changes to N-th (last) diagonal element
  // Is outsourced, since only the diagonal element is modified.
  if (*tbuff != 0.0 || wkvec->data[size - 1] != 0.0) {
    gsl_blas_drotg(tbuff, wkvec->data + (size - 1), cvec->data + i_N, svec->data + i_N);
    if ((temp = *tbuff) < 0.0) {
      *tbuff = -temp;
      cvec->data[i_N] = -cvec->data[i_N];
      svec->data[i_N] = -svec->data[i_N];
    } else if (temp == 0.0) {
      throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
    }
  } else {
    throw sgpp::base::data_exception("choleskyUpdate::Matrix not numerical positive definite");
  }
  gsl_vector_free(wkvec);
  gsl_vector_free(svec);
  gsl_vector_free(cvec);
#else
  throw base::not_implemented_exception("built withot GSL");
#endif /* USE_GSL */
}

// Implement cholesky Downdate for given Decomposition and update vector
void DBMatDMSChol::choleskyDowndate(sgpp::base::DataMatrix& decompMatrix,
                                    const sgpp::base::DataVector& downdate, bool do_cv) const {
#ifdef USE_GSL
  size_t size = decompMatrix.getNrows();

  if (downdate.getSize() != size) {
    throw sgpp::base::data_exception(
        "choleskyDowndate::Size of DecomposedMatrix and updateVector don´t "
        "match...");
  }
  // Create GSL matrix view for update procedures
  gsl_matrix_view m = gsl_matrix_view_array(decompMatrix.getPointer(), size, size);
  gsl_vector_const_view vvec =
      gsl_vector_const_view_array(downdate.getPointer(), downdate.getSize());

  // Define and declare Workingvector, Cosine- and Sinevector
  gsl_vector* wkvec = gsl_vector_calloc(downdate.getSize());
  gsl_vector* svec = gsl_vector_calloc(downdate.getSize());
  gsl_vector* cvec = gsl_vector_calloc(downdate.getSize());

  // Compute p (if not given)
  // Copy Values of updateVector into WorkingVector
  gsl_blas_dcopy(&vvec.vector, wkvec);

  // Solve La = downdate and save a in vvec.vector
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m.matrix, wkvec);

  // Generate Givens rotations
  double rho;
  gsl_blas_ddot(wkvec, wkvec, &rho);
  rho = 1 - rho;

  // Represents first index of downdate vector with entrie unequal to zero
  size_t cv_first_zero = 0;
  if (rho <= 0.0) {
    throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite");
  } else {
    rho = sqrt(rho);
    for (int i = static_cast<int>(size) - 1; i >= 0; i--) {
      // If this method is applied in cross-validation do_cv == true and
      // shortcuts can be used
      if (wkvec->data[i] == 0 && do_cv == true) {
        cv_first_zero = i + 1;
        break;
      }
      // Determine Givens rotation
      gsl_blas_drotg(&rho, wkvec->data + i, cvec->data + i, svec->data + i);
      // rho must remain positive
      if (rho < 0.0) {
        rho = -rho;
        cvec->data[i] = -cvec->data[i];
        svec->data[i] = -svec->data[i];
      }
    }
  }

  // rho should be 1 now
  gsl_vector_set_zero(wkvec);

  double* tbuff = m.matrix.data + ((size - 1) * (size + 1));

  // Apply calculated Givens rotations to current Cholesky factor
  for (size_t i = size - 1; i >= cv_first_zero; i--) {
    if (*tbuff <= 0.0) {
      throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite");
    }
    // Access corresponding columns
    gsl_vector_view mat_sub = gsl_matrix_subcolumn(&m.matrix, i, i, size - i);
    gsl_vector_view wkvec_sub = gsl_vector_subvector(wkvec, i, size - i);

    // Apply Givens rotation to mat_sub and wkvec_sub
    gsl_blas_drot(&wkvec_sub.vector, &mat_sub.vector, cvec->data[i], svec->data[i]);
    double diag = gsl_matrix_get(&m.matrix, i, i);
    // Ensure diagonal stays positive
    if (diag < 0.0) {
      rho = -1.0;
      gsl_matrix_set(&m.matrix, i, i, rho * diag);
    } else if (diag == 0.0) {
      throw sgpp::base::data_exception("choleskyDowndate::Matrix not numerical positive definite");
    }
    tbuff -= (size + 1);
  }

  // Workingvector should equal v now

  gsl_vector_free(wkvec);
  gsl_vector_free(svec);
  gsl_vector_free(cvec);
#else
  throw base::not_implemented_exception("built withot GSL");
#endif /* USE_GSL */
}

void DBMatDMSChol::choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix,
                                        double lambda_up) const {
  size_t size = decompMatrix.getNcols();

  sgpp::base::DataVector lambdaModification(size, 0.0);
  // In case lambda is increased apply Cholesky rank one updates
  if (lambda_up > 0) {
    for (size_t i = 0; i < size; i++) {
      lambdaModification.set(i, sqrt(lambda_up));
      choleskyUpdate(decompMatrix, lambdaModification, true);
      lambdaModification.setAll(0.0);
    }
  } else if (lambda_up < 0) {
    // In case lambda is decreased apply Cholesky rank one downdates
    for (size_t i = 0; i < size; i++) {
      lambdaModification.set(i, sqrt(fabs(lambda_up)));
      choleskyDowndate(decompMatrix, lambdaModification, true);
      lambdaModification.setAll(0.0);
    }
  }
}

void DBMatDMSChol::choleskyBackwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                         const sgpp::base::DataVector& y,
                                         sgpp::base::DataVector& alpha) const {
  size_t size = decompMatrix.getNcols();
  for (int i = static_cast<int>(size) - 1; i >= 0; i--) {
    alpha[i] = y[i];
    for (size_t j = i + 1; j < size; j++) {
      alpha[i] -= decompMatrix.get(j, i) * alpha[j];
    }
    alpha[i] /= decompMatrix.get(i, i);
  }
}

void DBMatDMSChol::choleskyForwardSolve(const sgpp::base::DataMatrix& decompMatrix,
                                        const sgpp::base::DataVector& b,
                                        sgpp::base::DataVector& y) const {
  size_t size = decompMatrix.getNcols();

  for (size_t i = 0; i < size; i++) {
    y[i] = b[i];
    for (size_t j = 0; j < i; j++) {
      y[i] -= decompMatrix.get(i, j) * y[j];
    }
    y[i] /= decompMatrix.get(i, i);
  }
}

}  // namespace datadriven
}  // namespace sgpp
