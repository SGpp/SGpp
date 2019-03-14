/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimationOnOffParallel.cpp
 *
 * Created on: Mar 13, 2019
 *     Author: Jan Schopohl
 */

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>

#include <list>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <string>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::RefinementFunctor;
using sgpp::base::RefinementFunctorType;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationOnOffParallel::ModelFittingDensityEstimationOnOffParallel(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation() {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

double ModelFittingDensityEstimationOnOffParallel::evaluate(const DataVector& sample) {
  // return online->eval(alpha, sample, *grid);
  return 0.0;
}

void ModelFittingDensityEstimationOnOffParallel::evaluate(DataMatrix& samples,
                                                          DataVector& results) {
  // online->eval(alpha, samples, results, *grid);
}

void ModelFittingDensityEstimationOnOffParallel::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fitParallel(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::fitParallel(DataMatrix& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
  auto& parallelConfig = this->config->getParallelConfig();

  // clear model
  reset();

  // build grid
  gridConfig.dim_ = newDataset.getNcols();
  std::cout << "Dataset dimension " << gridConfig.dim_ << std::endl;
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    // build surplus vector
    alpha = DataVector{grid->getSize()};
  }

  // Build the offline instance first
  DBMatOffline* offline = nullptr;

  // only run offline phase on master node
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    // Intialize database if it is provided
    if (!databaseConfig.filepath.empty()) {
      datadriven::DBMatDatabase database(databaseConfig.filepath);
      // Check if database holds a fitting lhs matrix decomposition
      if (database.hasDataMatrix(gridConfig, refinementConfig, regularizationConfig,
                                 densityEstimationConfig)) {
        std::string offlineFilepath = database.getDataMatrix(
            gridConfig, refinementConfig, regularizationConfig, densityEstimationConfig);
        offline = DBMatOfflineFactory::buildFromFile(offlineFilepath);
      }
    }

    // Build and decompose offline object if not loaded from database
    if (offline == nullptr) {
      // Build offline object by factory, build matrix and decompose
      offline = DBMatOfflineFactory::buildOfflineObject(
          gridConfig, refinementConfig, regularizationConfig, densityEstimationConfig);
      offline->buildMatrix(grid.get(), regularizationConfig);
      offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    }
  }

  // distribute offline, data and alpha
  /*DataMatrixDistributed dataDistributed(data.data(), processGrid, data.getNrows(),
  data.getNcols(), parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  DataVectorDistributed alphaDistributed(alpha.data(), processGrid, alpha.getSize(),
                                         parallelConfig.rowBlockSize_);*/

  // online phase
  online = std::unique_ptr<DBMatOnlineDE>{
      DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline, *grid, regularizationConfig.lambda_)};

  /*online->computeDensityFunctionParallel(alpha, newDataset, *grid,
                                         this->config->getDensityEstimationConfig(), true,
                                         this->config->getCrossvalidationConfig().enable_);*/
  online->setBeta(this->config->getLearnerConfig().beta);
  // online->normalize(alpha, *grid);
}

bool ModelFittingDensityEstimationOnOffParallel::refine(size_t newNoPoints,
                                                        std::list<size_t>* deletedGridPoints) {
  /*
  // Coarsening, remove idx from alpha
  if (deletedGridPoints != nullptr && deletedGridPoints->size() > 0) {
    // Restructure alpha
    std::vector<size_t> idxToDelete{std::begin(*deletedGridPoints), std::end(*deletedGridPoints)};
    alpha.remove(idxToDelete);
  }
  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // Refinement, expand alpha
  if (newNoPoints > oldNoPoints) {
    alpha.resizeZero(newNoPoints);
  }

  // Update online object: lhs, rhs and recompute the density function based on the b stored
  online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(), *grid,
                                          newNoPoints - oldNoPoints, *deletedGridPoints,
                                          config->getRegularizationConfig().lambda_);
  online->updateRhs(newNoPoints, deletedGridPoints);*/
  return true;
}

void ModelFittingDensityEstimationOnOffParallel::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::update(DataMatrix& newDataset) {
  /*if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Update the fit (streaming)
    online->computeDensityFunction(alpha, newDataset, *grid,
                                   this->config->getDensityEstimationConfig(), true,
                                   this->config->getCrossvalidationConfig().enable_);
    // online->normalize(alpha, *grid);
  }*/
}

bool ModelFittingDensityEstimationOnOffParallel::isRefinable() {
  if (grid != nullptr) {
    return online->getOfflineObject().isRefineable();
  }
  return false;
}

void ModelFittingDensityEstimationOnOffParallel::reset() {
  grid.reset();
  online.reset();
  refinementsPerformed = 0;
}

}  // namespace datadriven
}  // namespace sgpp
