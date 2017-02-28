// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>

namespace sgpp {
namespace datadriven {

DBMatOnline::DBMatOnline() : offlineObject(nullptr) {}

DBMatOnline::DBMatOnline(DBMatOffline* o) { readOffline(o); }

void DBMatOnline::setLambda(double lambda) {
  switch (offlineObject->getConfig()->decomp_type_) {
    case DBMatDecompostionType::DBMatDecompEigen:
    case DBMatDecompostionType::DBMatDecompChol:
      offlineObject->getConfig()->lambda_ = lambda;
      break;
    case DBMatDecompostionType::DBMatDecompIChol:
      sgpp::base::application_exception(
          "Lambda can not be changed for IChol - not implemented yet");
      break;
    case DBMatDecompostionType::DBMatDecompLU:
    default:
      throw sgpp::base::application_exception(
          "Lambda can not be changed in the online step for this decomposition "
          "type!");
  }
}

void DBMatOnline::readOffline(DBMatOffline* o) { offlineObject = o; }

DBMatOffline* DBMatOnline::getOffline() { return offlineObject; }

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
