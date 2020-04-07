// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <list>
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
      processGrid(std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                     config.getParallelConfig().processCols_)),
      alphaDistributed(processGrid, 1, 1) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

ModelFittingDensityEstimationOnOffParallel::ModelFittingDensityEstimationOnOffParallel(
    const FitterConfigurationDensityEstimation& config,
    std::shared_ptr<BlacsProcessGrid> processGrid)
    : ModelFittingDensityEstimation(),
      processGrid(processGrid),
      alphaDistributed(processGrid, 1, 1) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

double ModelFittingDensityEstimationOnOffParallel::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

void ModelFittingDensityEstimationOnOffParallel::evaluate(DataMatrix& samples,
                                                          DataVector& results) {
  auto& parallelConfig = this->config->getParallelConfig();
  DataVectorDistributed resultsDistributed(processGrid, results.size(),
                                           parallelConfig.rowBlockSize_);

  if (processGrid->isProcessInGrid()) {
    online->evalParallel(alpha, samples, resultsDistributed, *grid);
  }

  // only the master needs the result, as it calculates the score
  resultsDistributed.toLocalDataVector(results);
}

void ModelFittingDensityEstimationOnOffParallel::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::fit(DataMatrix& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& geometryConfig = this->config->getGeometryConfig();
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
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig, geometryConfig)};

  alpha = DataVector(grid->getSize());

  // Build the offline instance first
  DBMatOffline* offline = nullptr;

  // Intialize database if it is provided
  if (!databaseConfig.filePath_.empty()) {
    datadriven::DBMatDatabase database(databaseConfig.filePath_);
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

    // add supported parallel version of offline decompositions here
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
      offline->decomposeMatrixParallel(regularizationConfig, densityEstimationConfig, processGrid,
                                       parallelConfig);
      // Note: do NOT compute the explicit inverse here for SMW_ decompositions, as regularization
      // needs to be done first
    } else {
      offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    }
  }

  alphaDistributed =
      DataVectorDistributed(processGrid, grid->getSize(), parallelConfig.rowBlockSize_);

  // online phase
  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(
      *offline, *grid, regularizationConfig.lambda_, 0, densityEstimationConfig.decomposition_)};
#ifdef USE_SCALAPACK
  online->syncDistributedDecomposition(processGrid, parallelConfig);
  // in case of SMW decomposition type, the inverse of the matrix needs to be
  // computed also
  if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
      densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
    offline->compute_inverse_parallel(processGrid, parallelConfig);
  }
  online->computeDensityFunctionParallel(alphaDistributed, newDataset, *grid,
                                         this->config->getDensityEstimationConfig(),
                                         this->config->getParallelConfig(), processGrid, true,
                                         this->config->getCrossvalidationConfig().enable_);
#endif /* USE_SCALAPACK */
  online->setBeta(this->config->getLearnerConfig().learningRate_);

  alpha = alphaDistributed.toLocalDataVectorBroadcast();

  if (densityEstimationConfig.normalize_) {
    online->normalize(alpha, *grid);
  }
}

bool ModelFittingDensityEstimationOnOffParallel::adapt(size_t newNoPoints,
                                                       std::vector<size_t>& deletedGridPoints) {
  // Coarsening, remove idx from alpha
  if (deletedGridPoints.size() > 0) {
    // Restructure alpha
    alpha.remove(deletedGridPoints);
  }

  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // update the distributed vector
  auto& parallelConfig = this->config->getParallelConfig();
  alphaDistributed =
      DataVectorDistributed(alpha.data(), processGrid, alpha.size(), parallelConfig.rowBlockSize_);

  // Refinement, expand alpha
  if (newNoPoints > oldNoPoints) {
    alphaDistributed.resize(newNoPoints);
    alpha.resizeZero(newNoPoints);
  }

  // in case of SMW decomposition type, refinement is distributed/parallelized also with special
  // handling, therefore has to be passed additional parameters of parallelconfig
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
  if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
      densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
#ifdef USE_SCALAPACK
    sgpp::datadriven::DBMatOnlineDE_SMW* online_SMW_pointer;
    online_SMW_pointer = static_cast<sgpp::datadriven::DBMatOnlineDE_SMW*>(&*online);
    online_SMW_pointer->updateSystemMatrixDecompositionParallel(
        config->getDensityEstimationConfig(), *grid, newNoPoints - oldNoPoints, deletedGridPoints,
        config->getRegularizationConfig().lambda_, processGrid, parallelConfig);
#endif      /* USE_SCALAPACK */
  } else {  // every other decomposition type than SMW
    // Update online object: lhs, rhs and recompute the density function based on the b stored
    online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(), *grid,
                                            newNoPoints - oldNoPoints, deletedGridPoints,
                                            config->getRegularizationConfig().lambda_);
  }
  online->updateRhs(newNoPoints, deletedGridPoints);

  online->syncDistributedDecomposition(processGrid, parallelConfig);
  return true;
}

void ModelFittingDensityEstimationOnOffParallel::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityEstimationOnOffParallel::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Update the fit (streaming)
    online->computeDensityFunctionParallel(alphaDistributed, newDataset, *grid,
                                           this->config->getDensityEstimationConfig(),
                                           this->config->getParallelConfig(), processGrid, true,
                                           this->config->getCrossvalidationConfig().enable_);

    alpha = alphaDistributed.toLocalDataVectorBroadcast();

    if (this->config->getDensityEstimationConfig().normalize_) {
      online->normalize(alpha, *grid);
    }
  }
}

double ModelFittingDensityEstimationOnOffParallel::computeResidual(
    DataMatrix& validationData) const {
  DataVectorDistributed bValidation = online->computeWeightedBFromBatchParallel(
      validationData, *grid, this->config->getDensityEstimationConfig(),
      this->config->getParallelConfig(), this->processGrid, true);

  DataMatrixDistributed rMatrix = online->getOfflineObject().getUnmodifiedRDistributed(
      this->processGrid, this->config->getParallelConfig());

#ifdef USE_SCALAPACK
  // R * alpha - b_val
  rMatrix.mult(alphaDistributed, bValidation, false, 1.0, -1.0);

  DataVector result = bValidation.toLocalDataVectorBroadcast();

  return result.l2Norm();
#else
  throw base::not_implemented_exception("built withot ScaLAPACK");
#endif /* USE_SCALAPACK */
}

void ModelFittingDensityEstimationOnOffParallel::updateRegularization(double lambda) {
  if (grid != nullptr) {
    auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
    auto& parallelConfig = this->config->getParallelConfig();

    this->online->getOfflineObject().updateRegularizationParallel(lambda, this->processGrid,
                                                                  parallelConfig);

    // in case of SMW decomposition type, the inverse of the matrix needs to be computed also
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
        densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
      online->getOfflineObject().compute_inverse_parallel(processGrid, parallelConfig);
    }
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

void ModelFittingDensityEstimationOnOffParallel::resetTraining() {
  if (grid != nullptr) {
    alpha = DataVector(grid->getSize());
    this->online->resetTraining();
  }
}

std::shared_ptr<BlacsProcessGrid> ModelFittingDensityEstimationOnOffParallel::getProcessGrid()
    const {
  return processGrid;
}

}  // namespace datadriven
}  // namespace sgpp
