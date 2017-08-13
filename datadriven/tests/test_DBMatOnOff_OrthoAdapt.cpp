/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt_test.cpp
 */

#define BOOST_TEST_DYN_LINK
// #ifdef USE_GSL
#include <gsl/gsl_blas.h>

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE(OrthoAdapt_tests)

BOOST_AUTO_TEST_CASE(offline_object) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;
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

BOOST_AUTO_TEST_CASE(online_object) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 2;
  config.grid_level_ = 3;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::OrthoAdapt;

  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_base(config);
  offline_base.buildMatrix();
  offline_base.decomposeMatrix();

  // create offline object like it should be after refined/coarsened Pts
  config.grid_level_++;
  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_refined(config);
  offline_refined.buildMatrix();

  // create online object based on offline_base
  auto online = std::unique_ptr<sgpp::datadriven::DBMatOnlineDE>{
      sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(offline_base)};

  // gather points to refine from bigger lhs_matrix
  size_t oldSize = offline_base.getGrid().getStorage().getSize();
  size_t newSize = offline_refined.getGrid().getStorage().getSize();
  size_t numberOfPoints = newSize - oldSize;
  sgpp::base::DataMatrix refinePts(newSize, numberOfPoints);

  // write points to refine into refinePts matrix
  for (size_t i = oldSize; i < newSize; i++) {
    sgpp::base::DataVector refPt(newSize);
    offline_refined.getDecomposedMatrix(true).getColumn(i, refPt);
    refinePts.setColumn(i, refPt);
    // also adding lambdas to "diagonals"
    refinePts.set(i - 1, i - oldSize - 1, refinePts.get(i - 1, i - oldSize - 1) + config.lambda_);
  }
/*
  // online object performs refinement now
  sgpp::datadriven::DBMatOnlineDEOrthoAdapt* thisChildPtr =
      static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*online);
  std::vector<size_t> coarsePts = {};
  thisChildPtr->adapt(refinePts, coarsePts); */
}

// #else
// throw sgpp::base::algorithm_exception("USE_GSL is not set to true");
// #endif /* USE_GSL */
BOOST_AUTO_TEST_SUITE_END()
