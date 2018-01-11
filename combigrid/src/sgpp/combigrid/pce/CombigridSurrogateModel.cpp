// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
namespace combigrid {

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
