// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <assert.h>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {
DBMatObjectStore::DBMatObjectStore() : hasDatabase(false) {}

DBMatObjectStore::DBMatObjectStore(const std::string& filePath)
    : dbFilePath(filePath), hasDatabase(true) {}

int DBMatObjectStore::getObjectContainerIndex(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    bool searchBase) {
  // Iterate over objects to find match
  for (size_t i = 0; i < this->objects.size(); i++) {
    // Check whether configuration of the container matches
    if (this->objects[i].configMatches(gridConfig, geometryConfig, adaptivityConfig,
                                       regularizationConfig, densityEstimationConfig, searchBase))
      return i;
  }
  // If no object matches, return -1
  return -1;
}

const DBMatObjectStore::ObjectContainer& DBMatObjectStore::getObjectContainer(int index) const {
  return this->objects.at(index);
}

const DBMatOffline* DBMatObjectStore::getObject(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // Search for suitable offline object
  int index = this->getObjectContainerIndex(gridConfig, geometryConfig, adaptivityConfig,
                                            regularizationConfig, densityEstimationConfig);
  // If no suitable object is found, return nullptr
  if (index == -1) {
    return nullptr;
  }
  // If suitable object is found, return pointer to the object
  else {
    return &this->getObjectContainer(index).getOfflineObject();
  }
}

const DBMatOfflinePermutable* DBMatObjectStore::getBaseObject(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    sgpp::base::GeneralGridConfiguration& baseGridConfig) {
  // Search for suitable base offline object
  int index = this->getObjectContainerIndex(gridConfig, geometryConfig, adaptivityConfig,
                                            regularizationConfig, densityEstimationConfig, true);
  // If no suitable base object is found, return nullptr
  if (index == -1) {
    return nullptr;
    // If suitable base object is found, return pointer to the object and base config
  } else {
    const DBMatObjectStore::ObjectContainer& result = this->getObjectContainer(index);
    baseGridConfig = result.getGridConfig();
    return (const DBMatOfflinePermutable*)&result.getOfflineObject();
  }
}

void DBMatObjectStore::putObject(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    const DBMatOffline* object) {
  // Add a new object container to stored containers
  this->objects.push_back(ObjectContainer{gridConfig, geometryConfig, adaptivityConfig,
                                          regularizationConfig, densityEstimationConfig,
                                          std::unique_ptr<const DBMatOffline>(object)});
}

DBMatObjectStore::ObjectContainer::ObjectContainer(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    std::unique_ptr<const DBMatOffline> offlineObject)
    : gridConfig(gridConfig),
      geometryConfig(geometryConfig),
      adaptivityConfig(adaptivityConfig),
      regularizationConfig(regularizationConfig),
      densityEstimationConfig(densityEstimationConfig),
      // Transfers ownership to container
      offlineObject(std::move(offlineObject)) {}

const DBMatOffline& DBMatObjectStore::ObjectContainer::getOfflineObject() const {
  return *(this->offlineObject);
}

const sgpp::base::GeneralGridConfiguration& DBMatObjectStore::ObjectContainer::getGridConfig()
    const {
  return this->gridConfig;
}

bool DBMatObjectStore::ObjectContainer::configMatches(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    bool searchBase) {
  bool configMatch = true;
  // If searching for a base object, it is checked whether the level vector of the potential base
  // object is a permutation of the desired object's level vector
  if (searchBase) {
    // Permutation and bloe-up approach is only applicable to component grids
    if (gridConfig.generalType_ != sgpp::base::GeneralGridType::ComponentGrid)
      throw sgpp::base::not_implemented_exception(
          "Base object can only be found for component grids.");
    // Remove elements equal to 1 from the desired level vector to check for permutation
    std::vector<size_t> levelVecWithoutOnes =
        PermutationUtil::deleteOnesFromLevelVec(gridConfig.levelVector_);
    // Ceck whether permutation from the stored offline object's level vector to the normalized
    // desired level vector exists
    configMatch =
        PermutationUtil::isPermutation(this->gridConfig.levelVector_, levelVecWithoutOnes);
  } else {
    // If an identical offline object is searced, simply check equality of grid configuration
    configMatch = configMatch && gridConfig.level_ == this->gridConfig.level_ &&
                  gridConfig.levelVector_ == this->gridConfig.levelVector_ &&
                  gridConfig.dim_ == this->gridConfig.dim_;
  }
  // All other cofiguration entries have to be equal
  return configMatch && gridConfig.boundaryLevel_ == this->gridConfig.boundaryLevel_ &&
         gridConfig.type_ == this->gridConfig.type_ && gridConfig.type_ == this->gridConfig.type_ &&
         // adaptivity config
         adaptivityConfig.errorBasedRefinement == this->adaptivityConfig.errorBasedRefinement &&
         adaptivityConfig.errorBufferSize == this->adaptivityConfig.errorBufferSize &&
         adaptivityConfig.errorConvergenceThreshold ==
             this->adaptivityConfig.errorConvergenceThreshold &&
         adaptivityConfig.errorMinInterval == this->adaptivityConfig.errorMinInterval &&
         adaptivityConfig.levelPenalize == this->adaptivityConfig.levelPenalize &&
         adaptivityConfig.maxLevelType_ == this->adaptivityConfig.maxLevelType_ &&
         adaptivityConfig.noPoints_ == this->adaptivityConfig.noPoints_ &&
         adaptivityConfig.numRefinements_ == this->adaptivityConfig.numRefinements_ &&
         adaptivityConfig.percent_ == this->adaptivityConfig.percent_ &&
         adaptivityConfig.precomputeEvaluations == this->adaptivityConfig.precomputeEvaluations &&
         adaptivityConfig.refinementFunctorType == this->adaptivityConfig.refinementFunctorType &&
         adaptivityConfig.refinementPeriod == this->adaptivityConfig.refinementPeriod &&
         adaptivityConfig.scalingCoefficients == this->adaptivityConfig.scalingCoefficients &&
         adaptivityConfig.threshold_ == this->adaptivityConfig.threshold_ &&
         // remaining regularization config
         regularizationConfig.exponentBase_ == this->regularizationConfig.exponentBase_ &&
         regularizationConfig.l1Ratio_ == this->regularizationConfig.l1Ratio_ &&
         regularizationConfig.type_ == this->regularizationConfig.type_ &&
         // density estimation config
         densityEstimationConfig.decomposition_ == this->densityEstimationConfig.decomposition_ &&
         densityEstimationConfig.iCholSweepsDecompose_ ==
             this->densityEstimationConfig.iCholSweepsDecompose_ &&
         densityEstimationConfig.iCholSweepsRefine_ ==
             this->densityEstimationConfig.iCholSweepsRefine_ &&
         densityEstimationConfig.iCholSweepsSolver_ ==
             this->densityEstimationConfig.iCholSweepsSolver_ &&
         densityEstimationConfig.iCholSweepsUpdateLambda_ ==
             this->densityEstimationConfig.iCholSweepsUpdateLambda_ &&
         densityEstimationConfig.normalize_ == this->densityEstimationConfig.normalize_ &&
         densityEstimationConfig.type_ == this->densityEstimationConfig.type_;
}
}  // namespace datadriven
}  // namespace sgpp
