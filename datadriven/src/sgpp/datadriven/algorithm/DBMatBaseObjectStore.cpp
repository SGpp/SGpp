#include <assert.h>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatBaseObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
DBMatBaseObjectStore::DBMatBaseObjectStore(
    sgpp::base::AdaptivityConfiguration adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig)
    : hasDatabase(false) {
  this->adaptivityConfig = adaptivityConfig;
  this->regularizationConfig = regularizationConfig;
  this->densityEstimationConfig = densityEstimationConfig;
}

DBMatBaseObjectStore::DBMatBaseObjectStore(
    const std::string& filePath, sgpp::base::AdaptivityConfiguration adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig)
    : dbFilePath(filePath), hasDatabase(true) {
  this->adaptivityConfig = adaptivityConfig;
  this->regularizationConfig = regularizationConfig;
  this->densityEstimationConfig = densityEstimationConfig;
}

DBMatOfflinePermutable* DBMatBaseObjectStore::getOfflineObject(
    sgpp::base::GeneralGridConfiguration& gridConfig) {
  // base object that will be transformed
  DBMatOfflinePermutable* baseObject = nullptr;
  // grid configuration of the base object
  sgpp::base::GeneralGridConfiguration baseGridConfig;
  // try to find matching object locally
  int baseObjectContainerIndex = this->getBaseObjectContainerIndex(gridConfig);
  // if a suitable object exists, use it
  if (baseObjectContainerIndex != -1) {
    std::cout << "Obtained Base Object from local store" << std::endl;
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
    std::cout << "Obtained Base Object from database" << std::endl;
    // contruct the database
    DBMatDatabase db(dbFilePath);
    if (db.hasBaseDataMatrix(gridConfig, this->adaptivityConfig, this->regularizationConfig,
                             this->densityEstimationConfig)) {
      sgpp::base::GeneralGridConfiguration dbGridConfig;
      std::string objectFile =
          db.getBaseDataMatrix(gridConfig, this->adaptivityConfig, this->regularizationConfig,
                               this->densityEstimationConfig, dbGridConfig);
      // build offline object
      baseObject = (DBMatOfflinePermutable*)DBMatOfflineFactory::buildFromFile(objectFile);
      // grid config with 1 elemts remove from level vector                    onst must be
      baseGridConfig = PermutationUtil::getNormalizedConfig(dbGridConfig);
      // permutate base object to match cleaned level vec
      baseObject->permutateDecomposition(dbGridConfig, baseGridConfig);
      // store base objects. Ownership gets transfered to the object's ObjectContainer
      this->putBaseObject(gridConfig, std::unique_ptr<DBMatOfflinePermutable>(baseObject));
    }
  }
  // if the object is not available, build it
  if (baseObject == nullptr) {
    std::cout << "Constructing base object ..." << std::endl;
    // remove 1 elements
    baseGridConfig = PermutationUtil::getNormalizedConfig(gridConfig);

    baseObject = (DBMatOfflinePermutable*)DBMatOfflineFactory::buildOfflineObject(
        baseGridConfig, this->adaptivityConfig, this->regularizationConfig,
        this->densityEstimationConfig);
    // build grid
    std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(baseGridConfig.dim_));
    grid->getGenerator().anisotropicFull(baseGridConfig.levelVector_);
    // build matrix
    baseObject->buildMatrix(grid.get(), this->regularizationConfig);
    // decompose matrix
    baseObject->decomposeMatrix(this->regularizationConfig, this->densityEstimationConfig);

    // store base objects. Ownership gets transfered to the object's ObjectContainer
    this->putBaseObject(baseGridConfig, std::unique_ptr<DBMatOfflinePermutable>(baseObject));

    std::cout << "Constructed base object ..." << std::endl;
  }
  // apply permutation and dimension blow-up
  baseObject->permutateDecomposition(baseGridConfig, gridConfig);
  return baseObject;
}

int DBMatBaseObjectStore::getBaseObjectContainerIndex(
    sgpp::base::GeneralGridConfiguration& gridConfig) {
  // Iterate over objects to find match
  for (int i = 0; i < this->objects.size(); i++) {
    if (this->objects[i].configMatches(gridConfig)) return i;
  }
  // If no object matches, return nullopt
  return -1;
}

const DBMatBaseObjectStore::ObjectContainer& DBMatBaseObjectStore::getBaseObjectContainer(
    int index) const {
  return this->objects.at(index);
}

void DBMatBaseObjectStore::putBaseObject(sgpp::base::GeneralGridConfiguration& gridConfig,
                                         std::unique_ptr<DBMatOfflinePermutable> object) {
  this->objects.push_back(ObjectContainer{gridConfig, std::move(object)});
}

DBMatBaseObjectStore::ObjectContainer::ObjectContainer(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    std::unique_ptr<const DBMatOfflinePermutable> offlineObject)
    : gridConfig(gridConfig), offlineObject(std::move(offlineObject)) {}

const DBMatOfflinePermutable& DBMatBaseObjectStore::ObjectContainer::getOfflineObject() const {
  return *(this->offlineObject);
}

const sgpp::base::GeneralGridConfiguration& DBMatBaseObjectStore::ObjectContainer::getGridConfig()
    const {
  return this->gridConfig;
}

bool DBMatBaseObjectStore::ObjectContainer::configMatches(
    const sgpp::base::GeneralGridConfiguration& gridConfig) {
  // store is only implemented for component grids
  if (gridConfig.generalType_ != sgpp::base::GeneralGridType::ComponentGrid)
    throw sgpp::base::not_implemented_exception("Only implemented for component grids.");
  std::vector<size_t> levelVecWithoutOnes = PermutationUtil::deleteOnesFromLevelVec(gridConfig.levelVector_);
  // check if a permutation exists
  return PermutationUtil::isPermutation(this->gridConfig.levelVector_, levelVecWithoutOnes);
  // remaining grid config
  /*gridConfig.boundaryLevel_ == this->gridConfig.boundaryLevel_ &&
    gridConfig.type_ == this->gridConfig.type_ &&
    // adaptivity config
    adaptivityConfig.errorBasedRefinement == this->adaptivityConfig.errorBasedRefinement &&
    adaptivityConfig.errorBufferSize == this->adaptivityConfig.errorBufferSize &&
    adaptivityConfig.errorConvergenceThreshold ==  this->adaptivityConfig.errorConvergenceThreshold
    && adaptivityConfig.errorMinInterval == this->adaptivityConfig.errorMinInterval &&
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
    densityEstimationConfig.type_ == this->densityEstimationConfig.type_;*/
}
}  // namespace datadriven
}  // namespace sgpp
