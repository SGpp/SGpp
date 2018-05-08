/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineEigen.cpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#ifdef USE_GSL
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::application_exception;
using sgpp::base::data_exception;
using sgpp::base::OperationMatrix;

DBMatOfflineEigen::DBMatOfflineEigen() {}

sgpp::datadriven::DBMatOfflineEigen::DBMatOfflineEigen(const std::string& fileName)
    : DBMatOffline{fileName} {}

DBMatOffline* DBMatOfflineEigen::clone() { return new DBMatOfflineEigen{*this}; }

bool DBMatOfflineEigen::isRefineable() { return false; }

void DBMatOfflineEigen::decomposeMatrix(RegularizationConfiguration& regularizationConfig,
    DensityEstimationConfiguration& densityEstimationConfig) {
  if (isConstructed) {
    if (isDecomposed) {
      // Already decomposed => Do nothing
      return;
    }
    size_t n = lhsMatrix.getNrows();

    gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                              n);  // Create GSL matrix view for decomposition

    auto q = std::unique_ptr<gsl_matrix>{gsl_matrix_alloc(n, n)};  // Stores the eigenvectors
    auto e = std::unique_ptr<gsl_vector>{gsl_vector_alloc(n)};     // Stores the eigenvalues

    gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(&m.matrix, e.get(), q.get(), ws);
    gsl_eigen_symmv_free(ws);

    // Create an (n+1)*n matrix to store eigenvalues and -vectors:
    lhsMatrix = DataMatrix(n + 1, n);

    for (size_t r = 0; r < n; r++) {
      for (size_t c = 0; c < n; c++) {
        lhsMatrix.set(r, c, gsl_matrix_get(q.get(), r, c));
      }
    }
    for (size_t c = 0; c < n; c++) {
      lhsMatrix.set(n, c, gsl_vector_get(e.get(), c));
    }

    isDecomposed = true;
  } else {
    throw data_exception("Matrix has to be constructed before it can be decomposed!");
  }
}

sgpp::datadriven::MatrixDecompositionType DBMatOfflineEigen::getDecompositionType() {
  return sgpp::datadriven::MatrixDecompositionType::Eigen;
}

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */
