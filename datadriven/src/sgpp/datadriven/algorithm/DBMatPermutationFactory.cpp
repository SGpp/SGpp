#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

namespace sgpp {
namespace datadriven {

DBMatPermutationFactory::DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store)
    : store(store), hasDataBase(false) {}

DBMatPermutationFactory::DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store,
                                                 const std::string& dbFilePath)
    : store(store), dbFilePath(dbFilePath), hasDataBase(true) {}

DBMatOfflinePermutable& DBMatPermutationFactory::getPermutatedObject(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // grid configuration of the base object
  sgpp::base::GeneralGridConfiguration baseGridConfig;
  // base object that will be transformed
  const DBMatOfflinePermutable* baseObject =
      this->store->getBaseObject(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                densityEstimationConfig, baseGridConfig);

  if (baseObject == nullptr) {
    // if no suitable object exists locally, and a db path is given, search the db and store the
    // base object if found
    if (hasDataBase) {
      // contruct the database
      DBMatDatabase db(dbFilePath);
      if (db.hasBaseDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
                               densityEstimationConfig)) {
        sgpp::base::GeneralGridConfiguration dbGridConfig;
        std::string objectFile =
            db.getBaseDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
                                 densityEstimationConfig, dbGridConfig);
        // build offline object
        DBMatOfflinePermutable* newBaseObject =
            (DBMatOfflinePermutable*)DBMatOfflineFactory::buildFromFile(objectFile);
        // grid config with 1 elemts remove from level vector                    onst must be
        baseGridConfig = PermutationUtil::getNormalizedConfig(dbGridConfig);
        // permutate base object to match cleaned level vec
        newBaseObject->permuteDecomposition(dbGridConfig, baseGridConfig);

        baseObject = newBaseObject;
        // store base objects. Ownership gets transfered to the object's ObjectContainer
        this->store->putObject(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                              densityEstimationConfig,
                              std::unique_ptr<const DBMatOfflinePermutable>(baseObject));
      }
    }
    // if the object is not available, build it
    else {
      // remove 1 elements
      baseGridConfig = PermutationUtil::getNormalizedConfig(gridConfig);

      // If normalized config has dimension 0, component grid has 1 level vector
      if (baseGridConfig.dim_ == 0) {
        baseGridConfig = gridConfig;
      }
      // Instanciate base offline object
      DBMatOfflinePermutable* newBaseObject =
          (DBMatOfflinePermutable*)DBMatOfflineFactory::buildOfflineObject(
              baseGridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);

      // build grid with geometry config
      std::unique_ptr<Grid> grid;
      sgpp::datadriven::GridFactory gridFactory;
      sgpp::datadriven::StencilType stencilType = geometryConfig.stencilType;
      std::vector<int64_t> dim = geometryConfig.dim;

      // a regular sparse grid is created, if no geometryConfig is defined,
      if (stencilType == sgpp::datadriven::StencilType::None) {
        // interaction with size 0
        std::vector<std::vector<size_t>> interactions = std::vector<std::vector<size_t>>();
        grid = std::unique_ptr<Grid>{gridFactory.createGrid(baseGridConfig, interactions)};
      } else {
        grid = std::unique_ptr<Grid>{
            gridFactory.createGrid(baseGridConfig, gridFactory.getInteractions(stencilType, dim))};
      }

      // build matrix
      newBaseObject->buildMatrix(grid.get(), regularizationConfig);
      // decompose matrix
      newBaseObject->decomposeMatrix(regularizationConfig, densityEstimationConfig);
      // store base objects, ownership gets transfered to the object's ObjectContainer
      this->store->putObject(baseGridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                            densityEstimationConfig,
                            std::unique_ptr<const DBMatOfflinePermutable>(baseObject));
    }
  }
  // copy return instance from base object
  DBMatOfflinePermutable* returnObject = dynamic_cast<DBMatOfflinePermutable*>(baseObject->clone());
  // apply permutation and dimension blow-up
  returnObject->permuteDecomposition(baseGridConfig, gridConfig);
  return *returnObject;
}
}  // namespace datadriven
}  // namespace sgpp