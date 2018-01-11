// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
namespace combigrid {

void CombigridSurrogateModelConfiguration::setCombigridOperation(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> pcombigridOperation) {
  combigridOperation = pcombigridOperation;
}

void CombigridSurrogateModelConfiguration::setCombigridMultiOperation(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> pcombigridMultiOperation) {
  combigridMultiOperation = pcombigridMultiOperation;
}

void CombigridSurrogateModelConfiguration::setCombigridTensorOperation(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> pcombigridTensorOperation) {
  combigridTensorOperation = pcombigridTensorOperation;
}

// ---------------------------------------------------------------------------------

CombigridSurrogateModel::CombigridSurrogateModel(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : config(config) {
  // initialize number of dimensions
  if (config.combigridOperation != nullptr) {
    // make sure that the number of dimensions match
    numDims = config.combigridOperation->numDims();
  } else if (config.combigridMultiOperation != nullptr) {
    numDims = config.combigridMultiOperation->numDims();
  } else if (config.combigridTensorOperation != nullptr) {
    numDims = config.combigridTensorOperation->numDims();
  } else {
    throw sgpp::base::application_exception(
        "CombigridSurrogateModel: no operation is set in surrogate model config");
  }
}

CombigridSurrogateModel::~CombigridSurrogateModel() {}

sgpp::combigrid::CombigridSurrogateModelConfiguration& CombigridSurrogateModel::getConfig() {
  return config;
}

} /* namespace combigrid */
} /* namespace sgpp */
