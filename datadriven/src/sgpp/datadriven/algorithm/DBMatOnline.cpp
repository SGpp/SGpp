// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <list>
#include <vector>

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

DBMatOnline::DBMatOnline(DBMatOffline& o) : offlineObject{o} {}

void DBMatOnline::setLambda(double lambda) {
  /**
  switch (offlineObject.getDensityEstimationConfig().decomposition_) {
    case MatrixDecompositionType::Eigen:
    case MatrixDecompositionType::Chol:
    case MatrixDecompositionType::DenseIchol:
    case MatrixDecompositionType::OrthoAdapt:
      offlineObject.getRegularizationConfig().lambda_ = lambda;
      break;
    case MatrixDecompositionType::LU:
    default:
      throw application_exception(
          "Lambda can not be changed in the online step for this decomposition "
          "type!");
  }
  */
}

DBMatOffline& DBMatOnline::getOfflineObject() {
  return const_cast<DBMatOffline&>(
      static_cast<const DBMatOnline&>(*this).getOfflineObject());
}

const DBMatOffline& DBMatOnline::getOfflineObject() const {
  return offlineObject;
}

std::vector<size_t> DBMatOnline::updateSystemMatrixDecomposition(
    DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
    size_t numAddedGridPoints, std::vector<size_t>& deletedGridPointIndices,
    double lambda) {
  if (!getOfflineObject().isRefineable()) {
    throw base::not_implemented_exception(
        "Attempted to update system matrix on decomposition that doesn't "
        "support it.");
  }
  throw base::application_exception(
      "Decomposition reports refineable but does not override "
      "updateSystemMatrixDecomposition()");
}

}  // namespace datadriven
}  // namespace sgpp
