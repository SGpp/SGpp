// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
// #ifdef USE_GSL
#include <gsl/gsl_blas.h>

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <iomanip>
#include <string>
#include <vector>

// print datamatrices for debugging
static void printMatrix(sgpp::base::DataMatrix a) {
  for (size_t i = 0; i < a.getNrows(); i++) {
    for (size_t j = 0; j < a.getNcols(); j++) {
      std::cout << std::setprecision(5) << std::fixed << a.get(i, j) << "  ";
    }
    std::cout << std::endl;
  }
}

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
  // std::cout << "Created Offline Object: \nMatrix Dimension = " << n << std::endl;
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
  config.grid_dim_ = 1;
  config.grid_level_ = 2;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::OrthoAdapt;

  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_base(config);
  offline_base.buildMatrix();
  sgpp::base::DataMatrix copy_small_lhs(offline_base.getDecomposedMatrix(true));
  offline_base.decomposeMatrix();

  // create offline object like it should be after refined Pts
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
    refinePts.setColumn(i - oldSize, refPt);
    // also adding lambdas to "diagonals"
    refinePts.set(i, i - oldSize, refinePts.get(i, i - oldSize) + config.lambda_);
  }

  // std::cout << "bigger grid lhsMatrix is \n";
  // printMatrix(offline_refined.getDecomposedMatrix(true));

  // std::cout << "\nand the refined points are: \n";
  // printMatrix(refinePts);

  /**
   * Testing refinement now
   */
  std::cout << "Testing refinement now ..." << std::endl;
  sgpp::datadriven::DBMatOnlineDEOrthoAdapt* thisChildPtr =
      static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*online);

  // refines!
  std::vector<size_t> coarsenPts = {6, 5, 4, 3};
  thisChildPtr->adapt(refinePts, true, coarsenPts);

  // create arbitrary right side b
  sgpp::base::DataVector b(newSize);
  double z = 1;
  for (size_t i = 0; i < b.getSize(); i++) {
    b.set(i, z);
    z = z * z * 0.25;
  }

  // create storage for alpha
  sgpp::base::DataVector alpha(newSize, 0.0);
  sgpp::base::DataVector test_result(newSize, 0.0);

  // solve the created system
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  solver->solve(offline_base.getTinv(), offline_base.getQ(), thisChildPtr->getB(), b, alpha);

  // now to test, it should hold that (lhs + lambda*I) * alpha = b
  // where the matrix is the one of the created big grid and alpha was
  // obtained from refining the smaller grid up to match the big one, ... so:
  for (size_t i = 0; i < newSize; i++) {
    double value = offline_refined.getDecomposedMatrix(true).get(i, i);
    offline_refined.getDecomposedMatrix(true).set(i, i, value + config.lambda_);
  }
  offline_refined.getDecomposedMatrix(true).mult(alpha, test_result);

  // test the results agains the original values of b
  for (size_t i = 0; i < alpha.getSize(); i++) {
    BOOST_CHECK_SMALL(test_result.get(i) - b.get(i), 1e-10);
  }

  /**
   * Testing coarsening now
   */
  std::cout << "Testing coarsening now ..." << std::endl;

  // now the algorithm can coarse the same points that were refined to test coarsening
  // refinePts must be sordet accordingly, because coarsen indices have to be in descendand
  // order, and must correspont to the index of refinePts
  // reminder: coarsePts = {6, 5, 4, 3}
  sgpp::base::DataMatrix coarsenMatrix(newSize, numberOfPoints);
  sgpp::base::DataVector puf(newSize, 0.0);

  // 3 --> 0
  refinePts.getColumn(3, puf);
  coarsenMatrix.setColumn(0, puf);
  // 2 --> 1
  refinePts.getColumn(2, puf);
  coarsenMatrix.setColumn(1, puf);
  // 1 --> 2
  refinePts.getColumn(1, puf);
  coarsenMatrix.setColumn(2, puf);
  // 0 --> 3
  refinePts.getColumn(0, puf);
  coarsenMatrix.setColumn(3, puf);

  // coarsening!
  thisChildPtr->adapt(coarsenMatrix, false, coarsenPts);

  // adapt vectors to the smaller dimension
  b.resize(oldSize);
  alpha.resize(oldSize);
  test_result.resize(oldSize);

  // solve the created system
  solver->solve(offline_base.getTinv(), offline_base.getQ(), thisChildPtr->getB(), b, alpha);

  // this time test agains smaller lhs matrix, (lhs + lambda*I) * alpha = b
  for (size_t i = 0; i < oldSize; i++) {
    double value = copy_small_lhs.get(i, i);
    copy_small_lhs.set(i, i, value + config.lambda_);
  }
  copy_small_lhs.mult(alpha, test_result);

  // test the results agains the original values of b
  for (size_t i = 0; i < alpha.getSize(); i++) {
    BOOST_CHECK_SMALL(test_result.get(i) - b.get(i), 1e-10);
  }
}

// #endif /* USE_GSL */
BOOST_AUTO_TEST_SUITE_END()
