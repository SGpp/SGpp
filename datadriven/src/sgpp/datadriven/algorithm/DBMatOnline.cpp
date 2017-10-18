// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

DBMatOnline::DBMatOnline(DBMatOffline& o) : offlineObject{o} {}

void DBMatOnline::setLambda(double lambda) {
  switch (offlineObject.getConfig().decomp_type_) {
    case DBMatDecompostionType::Eigen:
    case DBMatDecompostionType::Chol:
    case DBMatDecompostionType::DenseIchol:
    case DBMatDecompostionType::OrthoAdapt:
      offlineObject.getConfig().lambda_ = lambda;
      break;
    case DBMatDecompostionType::LU:
    default:
      throw application_exception(
          "Lambda can not be changed in the online step for this decomposition "
          "type!");
  }
}

DBMatOffline& DBMatOnline::getOfflineObject() {
  return const_cast<DBMatOffline&>(static_cast<const DBMatOnline&>(*this).getOfflineObject());
}

const DBMatOffline& DBMatOnline::getOfflineObject() const { return offlineObject; }

void DBMatOnline::updateSystemMatrixDecomposition(size_t numAddedGridPoints,
                                     std::list<size_t> deletedGridPointIndices,
                                     double lambda) {
  if(!getOfflineObject().isRefineable()) {
    throw base::not_implemented_exception("Attempted to update system matrix on decomposition "
                                                  "that doesn't support it.");
  }
  throw base::application_exception("Decomposition reports refineable but does not "
                                            "override updateSystemMatrixDecomposition()");
}


}  // namespace datadriven
}  // namespace sgpp
