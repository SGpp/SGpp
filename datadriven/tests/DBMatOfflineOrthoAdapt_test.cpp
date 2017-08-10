/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineOrthoAdapt_test.cpp
 *
 *  Created on: 01.08.2017
 *  Author: Dmitrij Boschko
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(OrthoAdapt_test)

BOOST_AUTO_TEST_CASE(offline_case) {
  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 3;
  config.grid_level_ = 3;  // grid lvl = 1 --> error
  config.grid_type_ = sgpp::base::GridType::Linear;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.lambda_ = 0.0001;

  sgpp::datadriven::DBMatOfflineOrthoAdapt off_object(config);

  size_t n = off_object.getDimA();
  std::cout << "Testing hessenberg_decomposition: \ndim = " << n << std::endl;
  BOOST_CHECK_EQUAL(n, 31);  // 31 points for dim = 3 and lvl = 3

  std::cout << ".";
  off_object.hessenberg_decomposition();
  std::cout << ".";

  // checks Q for orthogonality
  sgpp::base::DataVector row(n, 0.0);
  sgpp::base::DataVector col(n, 0.0);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      off_object.getQ().getColumn(i, col);
      off_object.getQ().getColumn(j, row);  // getColumn equals getRow of Q_transpose
      const double value = row.dotProduct(col);
      if (i == j) {
        BOOST_CHECK_CLOSE(value, 1.0, 1e-2);
      } else {
        BOOST_CHECK_SMALL(value, 1e-2);
      }
    }
  }
  std::cout << ".passed!\n";

  // copy diag and subdiag for validation
  // todo: copy constructor doesn't work, --> manual copy
  sgpp::base::DataVector diag(n);
  sgpp::base::DataVector subdiag(n);

  for (size_t i = 0; i < n; i++) {
    diag.set(i, off_object.getDiag().get(i));
  }
  for (size_t i = 0; i < n - 1; i++) {
    subdiag.set(i, off_object.getSubDiag().get(i));
  }

  std::cout << "Testing invert_tridiag: \n.";
  off_object.invert_tridiag();
  std::cout << ".";

  // checks T_inv * T = Id
  for (size_t i = 0; i < n; i++) {
    // make row of T_inv
    off_object.getTinv().getRow(i, row);
    for (size_t j = 0; j < n; j++) {
      // std::cout << " \n\ni=" << i << "   j=" << j << std::endl;
      // make column of T
      col.setAll(0.0);
      if (j == 0) {
        col.set(j + 1, subdiag.get(j));
        col.set(j, diag.get(j));
      } else if (j == n - 1) {
        col.set(j - 1, subdiag.get(j));
        col.set(j, diag.get(j));
      } else {
        col.set(j - 1, subdiag.get(j - 1));
        col.set(j, diag.get(j));
        col.set(j + 1, subdiag.get(j));
      }
      const double value = row.dotProduct(col);
      if (i == j) {
        BOOST_CHECK_CLOSE(value, 1.0, 1e-2);
      } else {
        BOOST_CHECK_SMALL(value, 1e-2);
      }
    }
  }
  std::cout << ".passed!\n";
}

BOOST_AUTO_TEST_SUITE_END()
