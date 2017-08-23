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
#include <list>
#include <string>
#include <vector>

static size_t static_dim = 1;
static int static_lvl = 2;

// print datamatrices for debugging
static void printMatrix(sgpp::base::DataMatrix a) {
  for (size_t i = 0; i < a.getNrows(); i++) {
    for (size_t j = 0; j < a.getNcols(); j++) {
      std::cout << std::setprecision(10) << std::fixed << a.get(i, j) << "  ";
    }
    std::cout << std::endl;
  }
}

static void calc_A_inv(sgpp::base::DataMatrix Q, sgpp::base::DataMatrix T_inv,
                       sgpp::base::DataMatrix A) {
  gsl_matrix_view q_view = gsl_matrix_view_array(Q.getPointer(), Q.getNrows(), Q.getNcols());
  gsl_matrix_view t_view =
      gsl_matrix_view_array(T_inv.getPointer(), T_inv.getNrows(), T_inv.getNcols());
  gsl_matrix_view a_view = gsl_matrix_view_array(A.getPointer(), A.getNrows(), A.getNcols());

  gsl_matrix* interim = gsl_matrix_alloc(Q.getNrows(), Q.getNcols());

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &q_view.matrix, &t_view.matrix, 0.0, interim);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, interim, &q_view.matrix, 0.0, &a_view.matrix);
  std::cout << "\ncalculated inverse of A: \n";
  printMatrix(A);
}

BOOST_AUTO_TEST_SUITE(OrthoAdapt_tests)

BOOST_AUTO_TEST_CASE(offline_object) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = static_dim;
  config.grid_level_ = static_lvl;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;

  sgpp::datadriven::DBMatOfflineOrthoAdapt off_object(config);
  off_object.buildMatrix();
  std::cout << "testing matrix size " << off_object.getDimA() << std::endl;

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
        BOOST_CHECK_SMALL(value - 1.0, 1e-12);
      } else {
        BOOST_CHECK_SMALL(value, 1e-12);
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
        BOOST_CHECK_SMALL(value - 1.0, 1e-12);
      } else {
        BOOST_CHECK_SMALL(value, 1e-12);
      }
    }
  }
  delete[] tt;
}

BOOST_AUTO_TEST_CASE(solver_test) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = static_dim;
  config.grid_level_ = static_lvl;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;

  sgpp::datadriven::DBMatOfflineOrthoAdapt off_object(config);
  off_object.buildMatrix();
  sgpp::base::DataMatrix lhs_copy(off_object.getLhsMatrix_ONLY_FOR_TESTING());

  size_t n = lhs_copy.getNcols();
  gsl_matrix_view lhs_view = gsl_matrix_view_array(lhs_copy.getPointer(), n, n);
  for (size_t i = 0; i < n; i++) {
    double val = gsl_matrix_get(&lhs_view.matrix, i, i) + config.lambda_;
    gsl_matrix_set(&lhs_view.matrix, i, i, val);
  }
  std::cout << "lhs_matrix + lambda*I: \n";
  printMatrix(lhs_copy);

  sgpp::base::DataVector b_array(n);
  for (size_t i = 0; i < n; i++) {
    b_array.set(i, i);
  }
  gsl_vector_view b_view = gsl_vector_view_array(b_array.getPointer(), n);
  sgpp::base::DataVector b(b_array);

  gsl_linalg_HH_svx(&lhs_view.matrix, &b_view.vector);

  sgpp::base::DataMatrix B_dummy(1, 1);
  sgpp::base::DataVector alpha(n);

  // solve the created system
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  off_object.decomposeMatrix();
  solver->solve(off_object.getTinv(), off_object.getQ(), B_dummy, b, alpha);

  // compare alpha values
  std::cout << "alpha from gsl_solver = \n";
  for (size_t i = 0; i < n; i++) {
    std::cout << std::setprecision(20) << b_view.vector.data[i] << "   ";
  }
  std::cout << std::endl;
  std::cout << "alpha from own solver = \n";
  for (size_t i = 0; i < n; i++) {
    std::cout << std::setprecision(20) << alpha.get(i) << "   ";
  }
  std::cout << std::endl;

  // check alphas
  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_SMALL(alpha.get(i) - b_view.vector.data[i], 1e-4);  // checkrange = lambda
  }
}

BOOST_AUTO_TEST_CASE(online_object) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = static_dim;
  config.grid_level_ = static_lvl;
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;
  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::OrthoAdapt;

  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_base(config);
  offline_base.buildMatrix();
  sgpp::base::DataMatrix copy_small_lhs(offline_base.getLhsMatrix_ONLY_FOR_TESTING());
  sgpp::base::DataMatrix copy_small_lhs_for_solve(offline_base.getLhsMatrix_ONLY_FOR_TESTING());

  offline_base.decomposeMatrix();

  // create offline object like it should be after refined Pts
  config.grid_level_++;
  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_refined(config);
  offline_refined.buildMatrix();

  // create online object based on offline_base
  auto online_parent = std::unique_ptr<sgpp::datadriven::DBMatOnlineDE>{
      sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(offline_base)};
  sgpp::datadriven::DBMatOnlineDEOrthoAdapt* online =
      static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*online_parent);

  // gather points to refine from bigger lhs_matrix
  size_t oldSize = offline_base.getGrid().getStorage().getSize();
  size_t newSize = offline_refined.getGrid().getStorage().getSize();
  size_t numberOfPoints = newSize - oldSize;
  std::vector<sgpp::base::DataVector> refine_points = {};

  // write points to refine into refinePts matrix
  for (size_t i = oldSize; i < newSize; i++) {
    sgpp::base::DataVector refPt(newSize);
    offline_refined.getLhsMatrix_ONLY_FOR_TESTING().getColumn(i, refPt);
    refPt.set(i, refPt.get(i) + config.lambda_);
    online->add_new_refine_point(refPt);
  }

  std::cout << "smaller grid lhsMatrix is \n";
  printMatrix(copy_small_lhs);

  std::cout << "bigger grid lhsMatrix is \n";
  printMatrix(offline_refined.getLhsMatrix_ONLY_FOR_TESTING());

  std::cout << "\nand the refined points are: \n";
  for (size_t j = 0; j < newSize; j++) {
    for (size_t i = 0; i < numberOfPoints; i++) {
      std::cout << std::setprecision(7) << std::fixed
                << (*online->getRefinedPointsPointer())[i].get(j) << "  ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // refines!
  std::list<size_t> coarsenPtsDummy = {};
  online->adapt(numberOfPoints, coarsenPtsDummy, config.lambda_);

  // create arbitrary right side b
  sgpp::base::DataVector b(newSize);
  double z = 0;
  for (size_t i = 0; i < b.getSize(); i++) {
    b.set(i, z);
    z++;
  }

  // create storage for alpha
  sgpp::base::DataVector alpha(newSize, 0.0);
  sgpp::base::DataVector test_result(newSize, 0.0);

  // solve the created system
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b, alpha);

  // now to test, it should hold that (lhs + lambda*I) * alpha = b
  // where the matrix is the one of the created big grid and alpha was
  // obtained from refining the smaller grid up to match the big one, ... so:
  for (size_t i = 0; i < newSize; i++) {
    double value = offline_refined.getLhsMatrix_ONLY_FOR_TESTING().get(i, i);
    offline_refined.getLhsMatrix_ONLY_FOR_TESTING().set(i, i, value + config.lambda_);
  }
  std::cout << "lambda adapted on diagonal big lhs: " << std::endl;
  printMatrix(offline_refined.getLhsMatrix_ONLY_FOR_TESTING());
  offline_refined.getLhsMatrix_ONLY_FOR_TESTING().mult(alpha, test_result);

  std::cout << "alpha is: " << std::endl;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    std::cout << std::fixed << alpha.get(i) << "  ";
  }
  std::cout << std::endl;

  std::cout << "(R_big + lambda*I) * alpha is: " << std::endl;
  for (size_t i = 0; i < test_result.getSize(); i++) {
    std::cout << std::fixed << test_result.get(i) << "  ";
  }
  std::cout << std::endl;

  // test the results against the original values of b
  for (size_t i = 0; i < alpha.getSize(); i++) {
    BOOST_CHECK_SMALL(test_result.get(i) - b.get(i), 1e-10);
  }

  std::cout << "Testing coarsening now ..." << std::endl;

  // get indices for coarsen points
  std::list<size_t> coarsen_points = {};
  for (size_t i = oldSize; i < newSize; i++) {
    coarsen_points.push_back(i);
  }

  // coarsening!
  online->adapt(0, coarsen_points, config.lambda_);

  // adapt vectors to the smaller dimension
  b.resize(oldSize);
  alpha.resize(oldSize);
  test_result.resize(oldSize);

  // solve the created system
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b, alpha);

  std::cout << "alpha is: " << std::endl;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    std::cout << std::fixed << alpha.get(i) << "  ";
  }
  std::cout << std::endl;

  // sgpp::base::DataMatrix A_inv(oldSize, oldSize);
  // calc_A_inv(offline_base.getQ(), offline_base.getTinv(), A_inv);
  // std::cout << "\n\nkonkretes A_inv von der kleinen matrix: \n";
  // printMatrix(A_inv);
  // A_inv.mult(b, alpha);
  // std::cout << "but alpha should be ( " << alpha.getSize() << "-dim Vector)" << std::endl;
  // for (size_t i = 0; alpha.getSize(); i++) {
  // std::cout << std::setprecision(7) << std::fixed << alpha.get(i) << "  ";
  // }

  // this time test agains smaller lhs matrix, (lhs + lambda*I) * alpha = b
  for (size_t i = 0; i < oldSize; i++) {
    double value = copy_small_lhs.get(i, i);
    copy_small_lhs.set(i, i, value + config.lambda_);
  }
  std::cout << "small_lhs ist: \n";
  printMatrix(copy_small_lhs);
  copy_small_lhs.mult(alpha, test_result);

  std::cout << "und (R_small + lambda*I) * alpha ist\n";
  for (size_t i = 0; i < test_result.getSize(); i++) {
    std::cout << std::fixed << test_result.get(i) << "  ";
  }
  std::cout << std::endl;

  // test the results agains the original values of b

  for (size_t i = 0; i < alpha.getSize(); i++) {
    BOOST_CHECK_SMALL(test_result.get(i) - b.get(i), 1e-10);
  }
  std::cout << "DAS WAR TEST: dim = " << config.grid_dim_
            << "\n              lvl:  " << config.grid_level_ - 1 << " ---> " << config.grid_level_
            << std::endl;
  std::cout << "              size: " << offline_base.getDimA() << " ---> "
            << offline_refined.getDimA() << std::endl;
}

// #endif /* USE_GSL */
BOOST_AUTO_TEST_SUITE_END()
