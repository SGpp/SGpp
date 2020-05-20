// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <list>
#include <string>
#include <set>
#include <vector>

BOOST_AUTO_TEST_SUITE(OrthoAdapt_tests)

BOOST_AUTO_TEST_CASE(offline_object) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.0001;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  sgpp::datadriven::DBMatOfflineOrthoAdapt off_object;
  off_object.buildMatrix(grid.get(), regularizationConfig);

  size_t n = off_object.getGridSize();
  std::cout << "Grid size is " << grid->getStorage().getSize() << std::endl;
  std::cout << "Matrix size " << n << std::endl;

  std::cout << "Testing hessenberg_decomposition...\n";

  // allocating sub-, super- and diagonal vectors of T
  sgpp::base::DataVector diag(n);
  sgpp::base::DataVector subdiag(n - 1);

  off_object.hessenberg_decomposition(diag, subdiag);

  gsl_matrix_view q_view = gsl_matrix_view_array(off_object.getQ().getPointer(), n, n);
  gsl_matrix* test_matrix = gsl_matrix_alloc(n, n);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &q_view.matrix, &q_view.matrix, 0.0, test_matrix);

  // checks Q for orthogonality: Q*Q^t == I
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i == j) {
        BOOST_CHECK_SMALL(test_matrix->data[n * i + j] - 1.0, 1e-12);
      } else {
        BOOST_CHECK_SMALL(test_matrix->data[n * i + j], 1e-12);
      }
    }
  }

  // creating explicit T for testing
  sgpp::base::DataMatrix T(n, n, 0.0);
  for (size_t i = 0; i < n; i++) {
    // adding lambda to diagonal
    diag.set(i, diag.get(i) + regularizationConfig.lambda_);
    T.set(i, i, diag.get(i));
  }
  for (size_t i = 0; i < n - 1; i++) {
    T.set(i + 1, i, subdiag.get(i));
    T.set(i, i + 1, subdiag.get(i));
  }

  std::cout << "Testing invert_symmetric_tridiag...\n";
  off_object.invert_symmetric_tridiag(diag, subdiag);

  gsl_matrix_view T_view = gsl_matrix_view_array(T.getPointer(), n, n);
  gsl_matrix_view T_inv_view = gsl_matrix_view_array(off_object.getTinv().getPointer(), n, n);

  // basically does T*T_inv
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &T_view.matrix, &T_inv_view.matrix, 0.0,
                 test_matrix);

  // checks if T_inv is inverse: T_inv * T == I
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i == j) {
        BOOST_CHECK_SMALL(test_matrix->data[n * i + j] - 1.0, 1e-12);
      } else {
        BOOST_CHECK_SMALL(test_matrix->data[n * i + j], 1e-12);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(solver_test) {
  // clang-format off
  double Q_array[] = {1, 0,
                      0, -1};
  sgpp::base::DataMatrix Q(Q_array, 2, 2);

  double T_array[] = {1, 0,
                      0, 2};
  sgpp::base::DataMatrix T(T_array, 2, 2);

  double b_array[] = {1, 2, 3, 4, 5};
  sgpp::base::DataVector b_cut(b_array, 2);
  sgpp::base::DataVector b(b_array, 4);

  double B_array[] = {1, 2, 0, 0,
                      2, 1, 0, 0,
                      0, 0, 3, -1,
                      0, 0, -1, 3};
  sgpp::base::DataMatrix B(B_array, 4, 4);
  // clang-format on

  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();

  // first case: no refinement, B = 0
  sgpp::base::DataVector first_alpha(2);
  sgpp::base::DataMatrix B_dummy(1, 1);
  solver->solve(T, Q, B_dummy, b_cut, first_alpha);
  // hard-coded calculated alpha values
  BOOST_CHECK_EQUAL(first_alpha.get(0), 1);
  BOOST_CHECK_EQUAL(first_alpha.get(1), 4);

  // second case: simulated "refinement", B != 0, so: alpha = Q * T^{-1} * Q^t * b + B * b
  sgpp::base::DataVector second_alpha(4);
  solver->solve(T, Q, B, b, second_alpha);
  // hard-coded calculated alpha values
  BOOST_CHECK_EQUAL(second_alpha.get(0), 6);
  BOOST_CHECK_EQUAL(second_alpha.get(1), 8);
  BOOST_CHECK_EQUAL(second_alpha.get(2), 5);
  BOOST_CHECK_EQUAL(second_alpha.get(3), 9);
}

/**
 * This test creates the offline/online setup, gets somewhat arbitrary points to refine
 * from a bigger lhs matrix, and then performs 2x refinement, after that 2x coarsening
 * in a different order. Values of solved alpha always are checked in between
 */
BOOST_AUTO_TEST_CASE(online_object) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 1;
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 0.0001;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};

  // creating offline objects
  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_base;
  offline_base.buildMatrix(grid.get(), regularizationConfig);  // creating lhs matrix
  sgpp::base::DataMatrix lhs(offline_base.getLhsMatrix_ONLY_FOR_TESTING());
  sgpp::base::DataMatrix lhs_copy_small(offline_base.getLhsMatrix_ONLY_FOR_TESTING());
  offline_base.decomposeMatrix(regularizationConfig, densityEstimationConfig);
  // calculating Q and T^{-1}

  auto online_parent = std::unique_ptr<sgpp::datadriven::DBMatOnlineDE>{
      sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(
          offline_base, *grid, regularizationConfig.lambda_, 0.0,
          densityEstimationConfig.decomposition_)};
  sgpp::datadriven::DBMatOnlineDEOrthoAdapt* online =
      static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*online_parent);

  // creating offline object of one bigger lvl as source for points to refine
  gridConfig.level_++;
  std::unique_ptr<sgpp::base::Grid> grid_source = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};
  sgpp::datadriven::DBMatOfflineOrthoAdapt offline_source;
  offline_source.buildMatrix(&(*grid_source), regularizationConfig);

  // calculate sizes of old and new matrices
  size_t oldSize = grid->getStorage().getSize();
  size_t newSize = grid_source->getStorage().getSize();
  size_t numberOfNewPoints = newSize - oldSize;  // always even, due to grid middle point

  //############################################################################
  //
  // first refinement test: refine half of the points
  //
  //############################################################################
  std::cout << "Testing refinement ..." << std::endl;

  // build expected already refined lhs matrix
  sgpp::base::DataMatrix lhs_copy_big(lhs_copy_small);
  lhs_copy_big.resizeQuadratic(oldSize + numberOfNewPoints / 2);

  // write first half of points to container and the expected refined lhs matrix
  for (size_t i = oldSize; i < oldSize + numberOfNewPoints / 2; i++) {
    sgpp::base::DataVector refPt(newSize);
    offline_source.getLhsMatrix_ONLY_FOR_TESTING().getColumn(i, refPt);
    refPt.set(i, refPt.get(i) + regularizationConfig.lambda_);
    online->add_new_refine_point(refPt);  // pushes point to container

    // fill the corresponding rows/columns of matrix with the new data
    for (size_t j = 0; j < lhs_copy_big.getNrows(); j++) {
      lhs_copy_big.set(i, j, refPt.get(j));
      lhs_copy_big.set(j, i, refPt.get(j));
    }
  }

  // perform the refining
  std::vector<size_t> dummy_coarsen_points = {};
  online->sherman_morrison_adapt(numberOfNewPoints / 2, true, dummy_coarsen_points);

  // creating space for alpha and values for b
  sgpp::base::DataVector alpha_half_refined(oldSize + numberOfNewPoints / 2);
  sgpp::base::DataVector b_half_refined(oldSize + numberOfNewPoints / 2);
  double z = 0.0;
  for (size_t i = 0; i < b_half_refined.getSize(); i++) {
    b_half_refined.set(i, z);
    z++;
  }

  // create ortho_adapt solver and solve the system for alpha
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b_half_refined,
                alpha_half_refined);

  // compute expected alpha with gsl QR solver
  gsl_vector_view gsl_alpha_half_refined_view =
      gsl_vector_view_array(alpha_half_refined.getPointer(), oldSize + numberOfNewPoints / 2);

  gsl_vector_view gsl_b_half_refined_view =
      gsl_vector_view_array(b_half_refined.getPointer(), oldSize + numberOfNewPoints / 2);

  sgpp::base::DataMatrix lhs_half_refined(lhs_copy_big);
  gsl_matrix_view gsl_lhs_half_refined_view =
      gsl_matrix_view_array(lhs_half_refined.getPointer(), oldSize + numberOfNewPoints / 2,
                            oldSize + numberOfNewPoints / 2);

  gsl_vector* tau = gsl_vector_alloc(oldSize + numberOfNewPoints / 2);
  gsl_linalg_QR_decomp(&gsl_lhs_half_refined_view.matrix, tau);
  gsl_linalg_QR_solve(&gsl_lhs_half_refined_view.matrix, tau, &gsl_b_half_refined_view.vector,
                      &gsl_alpha_half_refined_view.vector);

  // check calculated alpha against expected alpha
  BOOST_CHECK_EQUAL(alpha_half_refined.getSize(), gsl_alpha_half_refined_view.vector.size);
  for (size_t i = 0; i < alpha_half_refined.getSize(); i++) {
    BOOST_CHECK_SMALL(alpha_half_refined.get(i) - gsl_alpha_half_refined_view.vector.data[i],
                      1e-10);
  }

  //############################################################################
  //
  // second refinement test: refine the other half of the points
  //
  //############################################################################

  // resize expected lhs matrix to full big size
  lhs_copy_big.resizeQuadratic(newSize);

  // write second half of points to container and the expected refined lhs matrix
  for (size_t i = oldSize + numberOfNewPoints / 2; i < newSize; i++) {
    sgpp::base::DataVector refPt(newSize);
    offline_source.getLhsMatrix_ONLY_FOR_TESTING().getColumn(i, refPt);
    refPt.set(i, refPt.get(i) + regularizationConfig.lambda_);
    online->add_new_refine_point(refPt);  // pushes point to container

    // fill the corresponding rows/columns of matrix with the new data
    for (size_t j = 0; j < lhs_copy_big.getNrows(); j++) {
      lhs_copy_big.set(i, j, refPt.get(j));
      lhs_copy_big.set(j, i, refPt.get(j));
    }
  }

  // perform the refining
  online->sherman_morrison_adapt(numberOfNewPoints / 2, true, dummy_coarsen_points);

  // creating space for alpha and values for b
  sgpp::base::DataVector alpha_full_refined(newSize);
  sgpp::base::DataVector b_full_refined(newSize);
  z = 0.0;
  for (size_t i = 0; i < b_full_refined.getSize(); i++) {
    b_full_refined.set(i, z);
    z++;
  }

  // solve the system for alpha
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b_full_refined,
                alpha_full_refined);

  // compute expected alpha with gsl QR solver
  gsl_vector_view gsl_alpha_full_refined_view =
      gsl_vector_view_array(alpha_full_refined.getPointer(), newSize);

  gsl_vector_view gsl_b_full_refined_view =
      gsl_vector_view_array(b_full_refined.getPointer(), newSize);

  sgpp::base::DataMatrix lhs_full_refined(lhs_copy_big);
  gsl_matrix_view gsl_lhs_full_refined_view =
      gsl_matrix_view_array(lhs_full_refined.getPointer(), newSize, newSize);

  tau = gsl_vector_alloc(newSize);
  gsl_linalg_QR_decomp(&gsl_lhs_full_refined_view.matrix, tau);
  gsl_linalg_QR_solve(&gsl_lhs_full_refined_view.matrix, tau, &gsl_b_full_refined_view.vector,
                      &gsl_alpha_full_refined_view.vector);

  // check calculated alpha against expected alpha
  BOOST_CHECK_EQUAL(alpha_full_refined.getSize(), gsl_alpha_full_refined_view.vector.size);
  for (size_t i = 0; i < alpha_full_refined.getSize(); i++) {
    BOOST_CHECK_SMALL(alpha_full_refined.get(i) - gsl_alpha_full_refined_view.vector.data[i],
                      1e-10);
  }

  //############################################################################
  //
  // first coarsening test: coarsen the first half of the points refined before
  //
  //############################################################################
  std::cout << "Testing coarsening ..." << std::endl;

  // create the list of indices which to coarsen
  std::vector<size_t> coarsen_indices_first_half;
  for (size_t i = oldSize; i < oldSize + numberOfNewPoints / 2; i++) {
    coarsen_indices_first_half.push_back(i);
  }

  // write second half of points to the expected refined lhs matrix
  // it should be the offline_base lhs matrix with added vectors from the second
  // batch of the refine points, because the first batch will be coarsened
  lhs_copy_big.resizeQuadratic(oldSize + numberOfNewPoints / 2);
  for (size_t i = oldSize + numberOfNewPoints / 2; i < newSize; i++) {
    sgpp::base::DataVector refPt(newSize);
    offline_source.getLhsMatrix_ONLY_FOR_TESTING().getColumn(i, refPt);
    refPt.set(i, refPt.get(i) + regularizationConfig.lambda_);
    // no adding to container here, since this is coarsening

    // erase the entries of the obtained rows/columns at coarsened indices and resize
    refPt.erase(refPt.begin() + oldSize, refPt.begin() + oldSize + numberOfNewPoints / 2);

    // fill the corresponding rows/columns of matrix with the new data
    for (size_t j = 0; j < lhs_copy_big.getNrows(); j++) {
      lhs_copy_big.set(i, j, refPt.get(j));
      lhs_copy_big.set(j, i, refPt.get(j));
    }
  }

  // perform the coarsening
  online->sherman_morrison_adapt(0, false, coarsen_indices_first_half);

  // solve the system for alpha
  // note: it doesn't matter that the entries of b didn't get erased on the correct positions
  //       because the validation works for any b, it just has to be the same
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b_half_refined,
                alpha_half_refined);

  // compute expected alpha with gsl QR solver
  sgpp::base::DataMatrix lhs_half_coarsened(lhs_copy_big);  // lhs_copy_big is coarsened matrix here
  gsl_matrix_view gsl_lhs_half_coarsened_view =
      gsl_matrix_view_array(lhs_half_coarsened.getPointer(), oldSize + numberOfNewPoints / 2,
                            oldSize + numberOfNewPoints / 2);

  tau = gsl_vector_alloc(oldSize + numberOfNewPoints / 2);
  gsl_linalg_QR_decomp(&gsl_lhs_half_coarsened_view.matrix, tau);
  gsl_linalg_QR_solve(&gsl_lhs_half_coarsened_view.matrix, tau, &gsl_b_half_refined_view.vector,
                      &gsl_alpha_half_refined_view.vector);

  // check calculated alpha against expected alpha
  BOOST_CHECK_EQUAL(alpha_full_refined.getSize(), gsl_alpha_full_refined_view.vector.size);
  for (size_t i = 0; i < alpha_full_refined.getSize(); i++) {
    BOOST_CHECK_SMALL(alpha_full_refined.get(i) - gsl_alpha_full_refined_view.vector.data[i],
                      1e-10);
  }

  //############################################################################
  //
  // second coarsening test: coarsen the rest of the points refined before
  //
  //############################################################################

  // create the list of indices which to coarsen
  std::vector<size_t> coarsen_indices_second_half;
  for (size_t i = oldSize + numberOfNewPoints / 2; i < newSize; i++) {
    coarsen_indices_second_half.push_back(i);
  }

  // in this case no need to create matrices, because the expected fully coarsened matrix
  // is the small lhs matrix of the initial offline grid, which was copied before

  // perform the coarsening
  online->sherman_morrison_adapt(0, false, coarsen_indices_first_half);

  // create alpha and b by resizing the old ones to the old size
  alpha_half_refined.resize(oldSize);
  b_half_refined.resize(oldSize);

  // solve the system for alpha, which was resized in the lines above
  solver->solve(offline_base.getTinv(), offline_base.getQ(), online->getB(), b_half_refined,
                alpha_half_refined);

  // compute expected alpha with gsl QR solver
  gsl_vector_view gsl_alpha_view = gsl_vector_view_array(alpha_half_refined.getPointer(), oldSize);

  gsl_vector_view gsl_b_view = gsl_vector_view_array(b_half_refined.getPointer(), oldSize);

  gsl_matrix_view gsl_lhs_view =
      gsl_matrix_view_array(lhs_copy_small.getPointer(), oldSize, oldSize);

  tau = gsl_vector_alloc(oldSize);
  gsl_linalg_QR_decomp(&gsl_lhs_view.matrix, tau);
  gsl_linalg_QR_solve(&gsl_lhs_view.matrix, tau, &gsl_b_view.vector, &gsl_alpha_view.vector);

  // check calculated alpha against expected alpha
  BOOST_CHECK_EQUAL(alpha_half_refined.getSize(), gsl_alpha_view.vector.size);
  for (size_t i = 0; i < alpha_half_refined.getSize(); i++) {
    BOOST_CHECK_SMALL(alpha_half_refined.get(i) - gsl_alpha_view.vector.data[i], 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
