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
    : ModelFittingDensityEstimation(),
      processGrid(std::make_shared<BlacsProcessGrid>()),
      alphaDistributed(processGrid, 1, 1) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

double ModelFittingDensityEstimationOnOffParallel::evaluate(const DataVector& sample) {
  throw sgpp::base::not_implemented_exception(
      "This function is not implemented in the parallelized version!");
  /*alpha = alphaDistributed.toLocalDataVector();
  double result = 0.0;

  // TODO(jan) only eval on one process (and broadcast) or an all processes?
  if (processGrid->getCurrentRow() == 0 && processGrid->getCurrentColumn() == 0) {
    result = online->eval(alpha, sample, *grid);

    // broadcast value to other processes so that every process returns the same value
    Cdgebs2d(processGrid->getContextHandle(), "All", "T", 1, 1, &result, 1);
  } else if (processGrid->isProcessInGrid()) {
    Cdgebr2d(processGrid->getContextHandle(), "All", "T", 1, 1, &result, 1, 0, 0);
  } else {
    std::cout << "Warning! Process not in the grid tried to call evaluate, invalid result will be "
                 "returned!"
              << std::endl;
  }
  return result;*/
}

void ModelFittingDensityEstimationOnOffParallel::evaluate(DataMatrix& samples,
                                                          DataVector& results) {
  alpha = alphaDistributed.toLocalDataVectorBroadcast();

  auto& parallelConfig = this->config->getParallelConfig();
  DataVectorDistributed resultsDistributed(processGrid, results.size(),
                                           parallelConfig.rowBlockSize_);

  if (processGrid->isProcessInGrid()) {
    online->evalParallel(alpha, samples, resultsDistributed, *grid);
  }

  resultsDistributed.toLocalDataVectorBroadcast(results);
}

void ModelFittingDensityEstimationOnOffParallel::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::fit(DataMatrix& newDataset) {
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

  // TODO(jan) only on master?
  // if (BlacsProcessGrid::getCurrentProcess() == 0) {
  // build surplus vector
  alpha = DataVector{grid->getSize()};
  //}

  // Build the offline instance first
  DBMatOffline* offline = nullptr;

  // TODO(jan) only run offline phase on master node??
  // if (BlacsProcessGrid::getCurrentProcess() == 0) {
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
  //}

  // distribute offline, data and alpha
  // TODO(jan) which processes get the data, distribute?
  /*DataMatrixDistributed dataDistributed(newDataset.data(), processGrid, newDataset.getNrows(),
                                        newDataset.getNcols(), parallelConfig.rowBlockSize_,
                                        parallelConfig.columnBlockSize_);*/

  alphaDistributed = DataVectorDistributed(alpha.data(), processGrid, grid->getSize(),
                                           parallelConfig.rowBlockSize_);

  // online phase
  online = std::unique_ptr<DBMatOnlineDE>{
      DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline, *grid, regularizationConfig.lambda_)};

  online->computeDensityFunctionParallel(
      alphaDistributed, newDataset, *grid, this->config->getDensityEstimationConfig(),
      this->config->getParallelConfig(), processGrid, true,  // TODO(jan) save_b true??
      this->config->getCrossvalidationConfig().enable_);
  online->setBeta(this->config->getLearnerConfig().beta);
  // online->normalize(alpha, *grid);
}

bool ModelFittingDensityEstimationOnOffParallel::refine(size_t newNoPoints,
                                                        std::list<size_t>* deletedGridPoints) {
  // alpha is currently updated on all processes TODO(jan) check performance
  alpha = alphaDistributed.toLocalDataVectorBroadcast();

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

  // update the distributed vector
  auto& parallelConfig = this->config->getParallelConfig();
  alphaDistributed =
      DataVectorDistributed(alpha.data(), processGrid, alpha.size(), parallelConfig.rowBlockSize_);

  // every process does exactly the same steps to update the matrix decomposition
  // TODO(jan) check performance

  // Update online object: lhs, rhs and recompute the density function based on the b stored
  online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(), *grid,
                                          newNoPoints - oldNoPoints, *deletedGridPoints,
                                          config->getRegularizationConfig().lambda_);
  online->updateRhs(newNoPoints, deletedGridPoints);
  return true;
}  // namespace datadriven

void ModelFittingDensityEstimationOnOffParallel::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // auto& parallelConfig = this->config->getParallelConfig();

    // distribute the new data TODO(jan) needed?
    /*DataMatrixDistributed dataDistributed(newDataset.data(), processGrid, newDataset.getNrows(),
                                          newDataset.getNcols(), parallelConfig.rowBlockSize_,
                                          parallelConfig.columnBlockSize_);*/

    // Update the fit (streaming)
    online->computeDensityFunctionParallel(
        alphaDistributed, newDataset, *grid, this->config->getDensityEstimationConfig(),
        this->config->getParallelConfig(), processGrid, true,  // TODO(jan) save_b true??
        this->config->getCrossvalidationConfig().enable_);
    // online->normalize(alpha, *grid);
  }
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
