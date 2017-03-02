// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

DBMatOnline::DBMatOnline(DBMatOffline& o) : offlineObject{o} {}

void DBMatOnline::setLambda(double lambda) {
  switch (offlineObject.getConfig().decomp_type_) {
    case DBMatDecompostionType::DBMatDecompEigen:
    case DBMatDecompostionType::DBMatDecompChol:
      offlineObject.getConfig().lambda_ = lambda;
      break;
    case DBMatDecompostionType::DBMatDecompIChol:
      application_exception("Lambda can not be changed for IChol - not implemented yet");
      break;
    case DBMatDecompostionType::DBMatDecompLU:
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

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
