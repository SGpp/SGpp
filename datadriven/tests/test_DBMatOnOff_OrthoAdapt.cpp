/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt_test.cpp
 */

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>

BOOST_AUTO_TEST_SUITE(OrthoAdapt_test)

BOOST_AUTO_TEST_CASE(offline_object) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;  // grid lvl = 1 --> error, grid lvl = 10 --> bad alloc
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;

  sgpp::datadriven::DBMatOfflineOrthoAdapt off_object(config);
  off_object.buildMatrix();

  size_t n = off_object.getDimA();
  std::cout << "Created Offline Object: \nMatrix Dimension = " << n << std::endl;
  std::cout << "Testing hessenberg_decomposition...\n";

  // allocating sub-, super- and diagonal vectors of T
  gsl_vector* gsl_diag = gsl_vector_alloc(n);
  gsl_vector* gsl_subdiag = gsl_vector_alloc(n - 1);

  off_object.hessenberg_decomposition(gsl_diag, gsl_subdiag);

  // checks Q for orthogonality: Q*Q^t == I
  sgpp::base::DataVector row(n, 0.0);
  sgpp::base::DataVector col(n, 0.0);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      off_object.getQ().getColumn(i, col);
      off_object.getQ().getColumn(j, row);  // getColumn equals getRow of Q_transpose
      const double value = row.dotProduct(col);
      if (i == j) {
        BOOST_CHECK_CLOSE(value, 1.0, 1e-8);
      } else {
        BOOST_CHECK_SMALL(value, 1e-10);
      }
    }
  }

  // creating explicit T for testing
  sgpp::base::DataMatrix T(n, n, 0.0);
  for (size_t i = 0; i < n; i++) {
    // adding lambda to diagonal
    gsl_vector_set(gsl_diag, i, gsl_vector_get(gsl_diag, i) + config.lambda_);
    T.set(i, i, gsl_vector_get(gsl_diag, i));
  }
  for (size_t i = 0; i < n - 1; i++) {
    T.set(i + 1, i, gsl_vector_get(gsl_subdiag, i));
    T.set(i, i + 1, gsl_vector_get(gsl_subdiag, i));
  }

  std::cout << "Testing invert_symmetric_tridiag...\n";
  off_object.invert_symmetric_tridiag(gsl_diag, gsl_subdiag);

  gsl_vector_free(gsl_diag);
  gsl_vector_free(gsl_subdiag);

  gsl_matrix_view T_view = gsl_matrix_view_array(T.getPointer(), n, n);
  gsl_matrix_view T_inv_view = gsl_matrix_view_array(off_object.getTinv().getPointer(), n, n);
  double* tt = new double[n * n];
  gsl_matrix_view tt_view = gsl_matrix_view_array(tt, n, n);

  // basically does T*T_inv
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &T_view.matrix, &T_inv_view.matrix, 0.0,
                 &tt_view.matrix);

  // checks if T_inv is inverse: T_inv * T == I
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      double value = tt[i * n + j];
      if (i == j) {
        BOOST_CHECK_CLOSE(value, 1.0, 1e-8);
      } else {
        BOOST_CHECK_SMALL(value, 1e-10);
      }
    }
  }
  delete[] tt;
}
#else
throw sgpp::base::algorithm_exception("USE_GSL is not set to true");
#endif /* USE_GSL */
BOOST_AUTO_TEST_SUITE_END()
