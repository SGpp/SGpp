/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimation.cpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

#include <string>
#include <list>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingBaseSingleGrid{}, refinementsPerformed{0}, offline{nullptr} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityEstimation::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityEstimation::evaluate(DataMatrix& samples, DataVector& results) {
    online->eval(alpha, samples, results, *grid);
}

void ModelFittingDensityEstimation::fit(Dataset& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();

  // clear model
  resetState();
  grid.reset();
  dataset = &newDataset;

  // build grid
  gridConfig.dim_ = dataset->getDimension();
  std::cout << "Dataset dimension " << gridConfig.dim_ << std::endl;
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  // grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  // build surplus vector
  alpha = DataVector{grid->getSize()};

  bool offlineInitialized = false;
  // Initialize database if provided (path is not empty)
  if (!databaseConfig.filepath.empty()) {
    datadriven::DBMatDatabase database(databaseConfig.filepath);
    // Check if database holds a fitting lhs matrix decomposition
    if (database.hasDataMatrix(gridConfig, refinementConfig, regularizationConfig,
      densityEstimationConfig)) {
      std::string offlineFilepath = database.getDataMatrix(gridConfig, refinementConfig,
         regularizationConfig, densityEstimationConfig);
      offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildFromFile(
         offlineFilepath)};
      offlineInitialized = true;
    }
  }

  // Build and decompose offline object if not loaded from database
  if (!offlineInitialized) {
    // Build offline object by factory, build matrix and decompose
    offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildOfflineObject(
       gridConfig, refinementConfig, regularizationConfig, densityEstimationConfig)};
    offline->buildMatrix(grid.get(), regularizationConfig);
    offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  }
  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline,
     *grid, regularizationConfig.lambda_)};

  online->computeDensityFunction(alpha, newDataset.getData(), *grid,
      this->config->getDensityEstimationConfig());
  online->normalize(alpha, *grid);
}

bool ModelFittingDensityEstimation::refine() {
  if (grid != nullptr && offline->isRefineable()) {
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      SurplusRefinementFunctor refinementFunctor(alpha, config->getRefinementConfig().noPoints_,
                                                 config->getRefinementConfig().threshold_);
      // refine grid
      auto oldNoPoints = grid->getSize();
      std::cout << "Old number points " << oldNoPoints << std::endl;
      grid->getGenerator().refine(refinementFunctor);
      auto newNoPoints = grid->getSize();
      std::cout << "New number points " << newNoPoints << std::endl;
      if (newNoPoints != oldNoPoints) {
        // Tell the SLE manager that the grid changed (for interal data structures)
        alpha.resizeZero(newNoPoints);
        std::list<size_t> deletedGridPoints;
        // TODO(roehner) enable coarsening
        online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(),
            *grid, newNoPoints - oldNoPoints, deletedGridPoints, online->getBestLambda());
        online->computeDensityFunction(alpha, dataset->getData(), *grid,
              this->config->getDensityEstimationConfig());
        online->normalize(alpha, *grid);
        refinementsPerformed++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }

  } else {
    throw application_exception(
        "ModelFittingDensityEstimation: Can't refine before initial grid is created");
    return false;
  }
  return false;
}

void ModelFittingDensityEstimation::update(Dataset& newDataset) {
  fit(newDataset);
}

void ModelFittingDensityEstimation::resetState() { refinementsPerformed = 0; }

}  // namespace datadriven
}  // namespace sgpp
