// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <string>

// #ifdef USE_GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace sgpp {
namespace datadriven {

class DBMatOfflineOrthoAdapt : public DBMatOffline {
 public:
  /**
   * Constructor
   * Builds DBMatOfflineOrthoAdapt Object from configuration
   *
   * @param config configuration for this offline object
   */
  explicit DBMatOfflineOrthoAdapt(const DBMatDensityConfiguration& config);

  /**
   * Constructor
   * Builds object from serialized offline object, stored in file
   *
   * @param fileName path to the file storing object
   */
  explicit DBMatOfflineOrthoAdapt(const std::string& fileName);

  DBMatOffline* clone();

  bool isRefineable();

  /**
   * Builds the right hand side matrix without the regularization term
   */
  void buildMatrix();

  /**
   * Decomposes and inverts the lhsMatrix of the offline object
   * (lhs + lambda*I)^-1 = Q * (T + lambda*I)^-1 * Q^t = Q * T_inv * Q^t
   *
   * The matrix lhsMatrix of the parent object will be altered during the process
   * uses: hessenberg_decomposition, invert_symmetric_tridiag
   */
  void decomposeMatrix();

  /**
   * Decomposes the lhsMatrix into lhs = Q * T * Q^t and stores the orthogonal
   * matrix Q into the member q_ortho_matrix. The information to reconstruct T
   * is written into diag and subdiag
   *
   * @param diag diagonal entries of T
   * @param subdiagonal and superdiagonal entries of T (symmetric)
   */
  void hessenberg_decomposition(gsl_vector* diag, gsl_vector* subdiag);

  /**
   * Inverts a symmetric tridiagonal matrix T, which is given in the form of
   * its diagonal and subdiagonal vectors. When finished, diag and subdiag no more
   * hold their initial values.
   *
   * @param diag diagonal entries of T
   * @param subdiag and superdiagonal entries of T (symmetric)
   */
  void invert_symmetric_tridiag(gsl_vector* diag, gsl_vector* subdiag);

  /**
   * Serialize the DBMatOfflineOrthoAdapt object
   * The lhsMatrix is stored in the form of compact tridiagonal decomposition,
   * which means the diagonal and subdiagonal of lhsMatrix are stored, and the
   * lower left part of the matrix holds the householder vectors.
   *
   * q_ortho_matrix and t_inv_tridiag
   * are also stored into the specified file, which is the explicit representation
   * of the decomposition needed for the online phase.
   *
   * @param fileName path where to store the file
   */
  void store(const std::string& fileName);

  // getter and setter
  size_t& getDimA() { return this->dim_a; };

  double getLambda() { return this->lambda; };

  sgpp::base::DataMatrix& getQ() { return this->q_ortho_matrix_; }

  sgpp::base::DataMatrix& getTinv() { return this->t_tridiag_inv_matrix_; }

 protected:
  size_t dim_a;                                  // quadratic matrix size of matrix to decompose
  double lambda;                                 // configuration
  sgpp::base::DataMatrix q_ortho_matrix_;        // orthogonal matrix
  sgpp::base::DataMatrix t_tridiag_inv_matrix_;  // inverse of the tridiag matrix
};
}  // namespace datadriven
}  // namespace sgpp
// #else
// throw sgpp::base::algorithm_exception("USE_GSL is not set to true");
// #endif /* USE_GSL */
