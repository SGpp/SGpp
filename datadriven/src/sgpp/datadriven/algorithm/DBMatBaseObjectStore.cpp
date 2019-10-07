#include <assert.h>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatBaseObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
DBMatBaseObjectStore::DBMatBaseObjectStore() : hasDatabase(false) {}

DBMatBaseObjectStore::DBMatBaseObjectStore(const std::string& filePath)
    : dbFilePath(filePath), hasDatabase(true) {}

DBMatOfflinePermutable* DBMatBaseObjectStore::getOfflineObject(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // base object that will be transformed
  DBMatOfflinePermutable* baseObject = nullptr;
  // grid configuration of the base object
  sgpp::base::GeneralGridConfiguration baseGridConfig;
  // try to find matching object locally
  int baseObjectContainerIndex = this->getBaseObjectContainerIndex(
      gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);
  // if a suitable object exists, use it
  if (baseObjectContainerIndex != -1) {
    // obtain the object container
    const DBMatBaseObjectStore::ObjectContainer& baseObjectContainer =
        this->getBaseObjectContainer(baseObjectContainerIndex);
    // copy the stored object in order to transform it
    baseObject = (DBMatOfflinePermutable*)baseObjectContainer.getOfflineObject().clone();
    // copy grid config from container
    baseGridConfig = baseObjectContainer.getGridConfig();
  }
  // if no suitable object exists locally, and a db path is given, search the db and store the base
  // object if found
  else if (hasDatabase) {
    // contruct the database
    DBMatDatabase db(dbFilePath);
    if (db.hasBaseDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
                             densityEstimationConfig)) {
      sgpp::base::GeneralGridConfiguration dbGridConfig;
      std::string objectFile =
          db.getBaseDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
                               densityEstimationConfig, dbGridConfig);
      // build offline object
      baseObject = (DBMatOfflinePermutable*)DBMatOfflineFactory::buildFromFile(objectFile);
      // grid config with 1 elemts remove from level vector                    onst must be
      baseGridConfig = PermutationUtil::getNormalizedConfig(dbGridConfig);
      // permutate base object to match cleaned level vec
      baseObject->permutateDecomposition(dbGridConfig, baseGridConfig);
      // store base objects. Ownership gets transfered to the object's ObjectContainer
      this->putBaseObject(gridConfig, adaptivityConfig, regularizationConfig,
                          densityEstimationConfig,
                          std::unique_ptr<DBMatOfflinePermutable>(baseObject));
    }
  }
  // if the object is not available, build it
  if (baseObject == nullptr) {
    // remove 1 elements
    baseGridConfig = PermutationUtil::getNormalizedConfig(gridConfig);
    baseObject = (DBMatOfflinePermutable*)DBMatOfflineFactory::buildOfflineObject(
        baseGridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);
    // build grid
    std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(baseGridConfig.dim_));
    // build matrix
    baseObject->buildMatrix(grid.get(), regularizationConfig);
    // decompose matrix
    baseObject->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    // store base objects. Ownership gets transfered to the object's ObjectContainer
    this->putBaseObject(gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig,
                        std::unique_ptr<DBMatOfflinePermutable>(baseObject));
  }
  // apply permutation and dimension blow-up
  baseObject->permutateDecomposition(baseGridConfig, gridConfig);
  return baseObject;
}

int DBMatBaseObjectStore::getBaseObjectContainerIndex(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // Iterate over objects to find match
  for (int i = 0; i < this->objects.size(); i++) {
    if (this->objects[i].configMatches(gridConfig, adaptivityConfig, regularizationConfig,
                                       densityEstimationConfig))
      return i;
  }
  // If no object matches, return nullopt
  return -1;
}

const DBMatBaseObjectStore::ObjectContainer& DBMatBaseObjectStore::getBaseObjectContainer(
    int index) const {
  return this->objects.at(index);
}

void DBMatBaseObjectStore::putBaseObject(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    std::unique_ptr<DBMatOfflinePermutable> object) {
  this->objects.push_back(ObjectContainer{gridConfig, adaptivityConfig, regularizationConfig,
                                          densityEstimationConfig, std::move(object)});
}

DBMatBaseObjectStore::ObjectContainer::ObjectContainer(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const RegularizationConfiguration& regularizationConfig,
    const DensityEstimationConfiguration& densityEstimationConfig,
    std::unique_ptr<const DBMatOfflinePermutable> offlineObject)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      regularizationConfig(regularizationConfig),
      densityEstimationConfig(densityEstimationConfig),
      offlineObject(std::move(offlineObject)) {}

const DBMatOfflinePermutable& DBMatBaseObjectStore::ObjectContainer::getOfflineObject() const {
  return *(this->offlineObject);
}

const sgpp::base::GeneralGridConfiguration& DBMatBaseObjectStore::ObjectContainer::getGridConfig()
    const {
  return this->gridConfig;
}

bool DBMatBaseObjectStore::ObjectContainer::configMatches(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // store is only implemented for component grids
  if (gridConfig.generalType_ != sgpp::base::GeneralGridType::ComponentGrid)
    throw sgpp::base::not_implemented_exception("Only implemented for component grids.");

  if (regularizationConfig.lambda_ != 0)
    throw sgpp::base::algorithm_exception(
        "Permutation is not possible with regularization applied. Lambda != 0.");
  // check if a permutation exists
  std::vector<size_t> levelVec(this->gridConfig.levelVector_);
  for (auto l : gridConfig.levelVector_) {
    // 1 elements are irrelevant
    if (l == 1) continue;
    // iterate through base level vector to find suitable element
    for (size_t i = 0; i < levelVec.size(); i++) {
      auto l_ = levelVec[i];
      // base objects are not allowed to have elements equal to 1.
      assert(l_ != 1);
      // if suitable elemet is found, remove it from the base vector by setting it to -1.
      if (l == l_) {
        levelVec[i] = -1;
        break;
      }
      // If no suitable element has been found, base level vec is no permutation
      else if (i == levelVec.size() - 1) {
        return false;
        break;
      }
    }
  }

  // base vec is permutation if all elements have been set to -1
  for (auto l : levelVec) {
    if (l != -1) return false;
  }

  // check if other configs match
  return  // remaining grid config
      gridConfig.boundaryLevel_ == this->gridConfig.boundaryLevel_ &&
      gridConfig.type_ == this->gridConfig.type_ &&
      // adaptivity config
      adaptivityConfig.errorBasedRefinement == this->adaptivityConfig.errorBasedRefinement &&
      adaptivityConfig.errorBufferSize == this->adaptivityConfig.errorConvergenceThreshold &&
      adaptivityConfig.levelPenalize == this->adaptivityConfig.levelPenalize &&
      adaptivityConfig.maxLevelType_ == this->adaptivityConfig.maxLevelType_ &&
      adaptivityConfig.noPoints_ == this->adaptivityConfig.numRefinements_ &&
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
