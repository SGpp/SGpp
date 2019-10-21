#include <assert.h>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
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
  for (int i = 0; i < this->objects.size(); i++) {
    if (this->objects[i].configMatches(gridConfig, searchBase)) return i;
  }
  // If no object matches, return nullopt
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
  int index = this->getObjectContainerIndex(gridConfig, geometryConfig, adaptivityConfig,
                                            regularizationConfig, densityEstimationConfig);

  if (index == -1)
    return nullptr;
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
  int index = this->getObjectContainerIndex(gridConfig, geometryConfig, adaptivityConfig,
                                            regularizationConfig, densityEstimationConfig, true);

  if (index == -1) {
    return nullptr;
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
    std::unique_ptr<const DBMatOffline> object) {
  this->objects.push_back(ObjectContainer{gridConfig, geometryConfig, adaptivityConfig,
                                          regularizationConfig, densityEstimationConfig,
                                          std::move(object)});
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
      offlineObject(std::move(offlineObject)) {}

const DBMatOffline& DBMatObjectStore::ObjectContainer::getOfflineObject() const {
  return *(this->offlineObject);
}

const sgpp::base::GeneralGridConfiguration& DBMatObjectStore::ObjectContainer::getGridConfig()
    const {
  return this->gridConfig;
}

bool DBMatObjectStore::ObjectContainer::configMatches(
    const sgpp::base::GeneralGridConfiguration& gridConfig, bool searchBase) {
  bool configMatch = true;
  if (searchBase) {
    if (gridConfig.generalType_ != sgpp::base::GeneralGridType::ComponentGrid)
      throw sgpp::base::not_implemented_exception(
          "Base object can only be found for component grids.");
    std::vector<size_t> levelVecWithoutOnes =
        PermutationUtil::deleteOnesFromLevelVec(gridConfig.levelVector_);
    // check if a permutation exists
    configMatch =
        PermutationUtil::isPermutation(this->gridConfig.levelVector_, levelVecWithoutOnes);
  } else {
    configMatch = configMatch && gridConfig.level_ == this->gridConfig.level_ &&
                  gridConfig.levelVector_ == this->gridConfig.levelVector_ &&
                  gridConfig.dim_ == this->gridConfig.dim_;
  }
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
