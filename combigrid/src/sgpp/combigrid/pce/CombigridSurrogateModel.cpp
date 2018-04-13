// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
namespace combigrid {

CombigridSurrogateModelConfiguration::CombigridSurrogateModelConfiguration()
    : type(CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION),
      pointHierarchies(0),
      storage(nullptr),
      levelManager(nullptr),
      levelStructure(nullptr),
      basisFunction(nullptr),
      basisFunctions(0),
      tensorOperation(nullptr),
      bounds(0),
      degree(3),
      coefficientStorage(nullptr),
      weightFunctions(0),
      enableLevelManagerStatsCollection(false),
      numDimensions(0) {}

CombigridSurrogateModelConfiguration::~CombigridSurrogateModelConfiguration() {}

void CombigridSurrogateModelConfiguration::loadFromCombigridOperation(
    std::shared_ptr<CombigridOperation> op, bool loadLevelStructure) {
  storage = op->getStorage();
  pointHierarchies = op->getPointHierarchies();
  numDimensions = pointHierarchies.size();
  if (loadLevelStructure) {
    levelStructure = op->getLevelManager()->getLevelStructure();
  }
}

void CombigridSurrogateModelConfiguration::loadFromCombigridOperation(
    std::shared_ptr<CombigridMultiOperation> op, bool loadLevelStructure) {
  storage = op->getStorage();
  pointHierarchies = op->getPointHierarchies();
  numDimensions = pointHierarchies.size();
  if (loadLevelStructure) {
    levelStructure = op->getLevelManager()->getLevelStructure();
  }
}

void CombigridSurrogateModelConfiguration::loadFromCombigridOperation(
    std::shared_ptr<CombigridTensorOperation> op, bool loadLevelStructure) {
  storage = op->getStorage();
  pointHierarchies = op->getPointHierarchies();
  tensorOperation = op;
  numDimensions = pointHierarchies.size();
  if (loadLevelStructure) {
    levelStructure = op->getLevelManager()->getLevelStructure();
  }
}

CombigridSurrogateModel::CombigridSurrogateModel(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : config(config), numDims(0) {
  // initialize number of dimensions
  if (config.pointHierarchies.size() > 0) {
    numDims = config.pointHierarchies.size();
  } else {
    throw sgpp::base::application_exception(
        "CombigridSurrogateModel: number of dimensions is unknown. Setting the point hierarchies "
        "is required.");
  }

  if ((config.basisFunctions.size() > 0 && config.basisFunctions.size() != numDims) ||
      (config.weightFunctions.size() > 0 && config.weightFunctions.size() != numDims) ||
      (config.bounds.size() > 0 && static_cast<size_t>(config.bounds.size() / 2) != numDims)) {
    throw sgpp::base::application_exception(
        "CombigridSurrogateModel: number of dimensions do not match.");
  }
}

CombigridSurrogateModel::~CombigridSurrogateModel() {}

sgpp::combigrid::CombigridSurrogateModelConfiguration& CombigridSurrogateModel::getConfig() {
  return config;
}

} /* namespace combigrid */
} /* namespace sgpp */
