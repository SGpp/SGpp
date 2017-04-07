/*
 * DBMatOfflineGE.cpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

#include <list>

namespace sgpp {
namespace datadriven {

using sgpp::base::operation_exception;
using sgpp::base::DataMatrix;

DBMatOfflineGE::DBMatOfflineGE(const DBMatDensityConfiguration& oc) : DBMatOffline(oc) {}

void DBMatOfflineGE::buildMatrix() {
  // build matrix
  DBMatOffline::buildMatrix();

  // then add regularization term
  auto size = grid->getStorage().getSize();

  // Construct matrix lambda * C (just use identity for C)
  DataMatrix lambdaC(size, size);
  if (config.regularization_ == RegularizationType::Identity) {
    lambdaC.setAll(0.);
    for (size_t i = 0; i < size; i++) {
      lambdaC.set(i, i, config.lambda_);
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
