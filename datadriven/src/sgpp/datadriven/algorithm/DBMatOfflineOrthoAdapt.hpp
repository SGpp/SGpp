/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt.hpp
 *
 *  Created on: 01.08.2017
 *  Author: Dmitrij Boschko
 */

// #ifdef USE_GSL

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace sgpp {
namespace datadriven {

/*
 * Implementation of offline phase, which uses orthogonal decomposition
 * (i.e. hessenberg decomposition) to obtain Q*T_inv*Q_transpose.
 * The decomposition lhs = Q*T*Q_t, and the inversion of T -> T_inv, are both
 * part of this class. --> (lhs + lambda*I)_inv = Q*T_inv*Q_transpose
 */
class DBMatOfflineOrthoAdapt : public DBMatOffline {
 public:
  explicit DBMatOfflineOrthoAdapt(const DBMatDensityConfiguration& config);

  explicit DBMatOfflineOrthoAdapt(const std::string& fileName);

  DBMatOffline* clone();

  bool isRefineable();

  void buildMatrix();

  void decomposeMatrix();

  /*
   * Decomposes an s.p.d. (symmetric and positive definite) matrix into the
   * factorization Q*T*Q_transpose, where Q is orthogonal and T is a symmetric
   * and tridiagonal matrix, which is stored into two vectors "diag" and "subdiag"
   */
  void hessenberg_decomposition();

  /*
   * Adds lambda*I to T and inverts it. This is equivalent with inverting
   * (lhs + lambda*I), because
   * (lhs + lambda*I)_inv = Q*(T + lambda*I)_inv*Q_transpose
   */
  void invert_tridiag();

  size_t& getDimA() { return this->dim_a; };
  void setDimA(size_t n) { this->dim_a = n; }

  sgpp::base::DataMatrix& getQ() { return this->q_ortho_matrix_; }

  sgpp::base::DataMatrix& getTinv() { return this->t_tridiag_inv_matrix_; }

  sgpp::base::DataVector& getDiag() { return this->diag_; }

  sgpp::base::DataVector& getSubDiag() { return this->subdiag_; }

  void setLHS(sgpp::base::DataMatrix& a) {
    for (size_t i = 0; i < dim_a; i++) {
      for (size_t j = 0; j < dim_a; j++) {
        this->lhsMatrix.set(i, j, a.get(i, j));
      }
    }
  }

 protected:
  size_t dim_a;
  sgpp::base::DataMatrix q_ortho_matrix_;
  sgpp::base::DataMatrix t_tridiag_inv_matrix_;
  sgpp::base::DataVector diag_;
  sgpp::base::DataVector subdiag_;
  double lambda;
};
}  // namespace datadriven
}  // namespace sgpp
// #endif /* USE_GSL */
