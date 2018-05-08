/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineGE.cpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#ifdef USE_GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#endif /* USE_GSL */

#include <list>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::operation_exception;
using sgpp::base::application_exception;
using sgpp::base::DataMatrix;

DBMatOfflineGE::DBMatOfflineGE() : DBMatOffline() {}


void DBMatOfflineGE::buildMatrix(Grid* grid, RegularizationConfiguration& regularizationConfig) {
  // build matrix
  DBMatOffline::buildMatrix(grid, regularizationConfig);

  // then add regularization term
  auto size = grid->getStorage().getSize();

  // Construct matrix lambda * C (just use identity for C)
  DataMatrix lambdaC(size, size);
  if (regularizationConfig.type_ == RegularizationType::Identity) {
    lambdaC.setAll(0.);
    for (size_t i = 0; i < size; i++) {
      lambdaC.set(i, i, regularizationConfig.lambda_);
    }
  } else {
    throw operation_exception("Unsupported regularization type");
  }

  // Compute A + lambda * C:
  lhsMatrix.add(lambdaC);

  isConstructed = true;
}

} /* namespace datadriven */
} /* namespace sgpp */
