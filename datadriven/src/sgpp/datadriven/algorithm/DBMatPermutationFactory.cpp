// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <set>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

DBMatPermutationFactory::DBMatPermutationFactory() : store(nullptr), hasDataBase(false) {}

DBMatPermutationFactory::DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store)
    : store(store), hasDataBase(false) {}

DBMatPermutationFactory::DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store,
                                                 const std::string& dbFilePath)
    : store(store), hasDataBase(true), dbFilePath(dbFilePath) {}

DBMatOfflinePermutable* DBMatPermutationFactory::getPermutedObject(
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
            dynamic_cast<DBMatOfflinePermutable*>(DBMatOfflineFactory::buildFromFile(objectFile));
        // grid config with 1 elemts remove from level vector                    onst must be
        baseGridConfig = PermutationUtil::getNormalizedConfig(dbGridConfig);
        // permutate base object to match cleaned level vec
        newBaseObject->permuteDecomposition(dbGridConfig, baseGridConfig);

        baseObject = newBaseObject;

        // store base objects. Ownership gets transfered to the object's ObjectContainer
        this->store->putObject(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                               densityEstimationConfig, baseObject);
      }
    } else {  // if the object is not available, build it
      // remove 1 elements
      baseGridConfig = PermutationUtil::getNormalizedConfig(gridConfig);

      // If normalized config has dimension 0, component grid has 1 level vector
      if (baseGridConfig.dim_ == 0) {
        baseGridConfig = gridConfig;
      }
      // Instanciate base offline object
      DBMatOfflinePermutable* newBaseObject =
          dynamic_cast<DBMatOfflinePermutable*>(DBMatOfflineFactory::buildOfflineObject(
              baseGridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig));

      // build grid with geometry config
      std::unique_ptr<Grid> grid;
      sgpp::datadriven::GridFactory gridFactory;

      // a regular sparse grid is created, if no geometryConfig is defined,
      if (geometryConfig.stencils_.empty()) {
        // interaction with size 0
        std::set<std::set<size_t>> interactions = std::set<std::set<size_t>>();
        grid = std::unique_ptr<Grid>{gridFactory.createGrid(baseGridConfig, interactions)};
      } else {
        grid = std::unique_ptr<Grid>{
            gridFactory.createGrid(baseGridConfig, gridFactory.getInteractions(geometryConfig))};
      }

      // build matrix
      newBaseObject->buildMatrix(grid.get(), regularizationConfig);
      // decompose matrix
      newBaseObject->decomposeMatrix(regularizationConfig, densityEstimationConfig);

      baseObject = newBaseObject;

      // store base objects, ownership gets transfered to the object's ObjectContainer
      this->store->putObject(baseGridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                             densityEstimationConfig, baseObject);
    }
  }
  // copy return instance from base object
  DBMatOfflinePermutable* returnObject = dynamic_cast<DBMatOfflinePermutable*>(baseObject->clone());
  // apply permutation and dimension blow-up
  returnObject->permuteDecomposition(baseGridConfig, gridConfig);
  return returnObject;
}
}  // namespace datadriven
}  // namespace sgpp
