// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOffParallel.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>
#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
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

ModelFittingDensityDerivativeEstimationOnOffParallel::
    ModelFittingDensityDerivativeEstimationOnOffParallel(
        const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation(),
      processGrid(std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                     config.getParallelConfig().processCols_)),
      alphaDistributed(processGrid, 1, 1) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  this->hasObjectStore = false;
}

ModelFittingDensityDerivativeEstimationOnOffParallel::
    ModelFittingDensityDerivativeEstimationOnOffParallel(
        const FitterConfigurationDensityEstimation& config,
        std::shared_ptr<DBMatObjectStore> objectStore,
        std::shared_ptr<BlacsProcessGrid> processGrid)
    : ModelFittingDensityEstimation(),
      objectStore(objectStore),
      hasObjectStore(false),
      processGrid(processGrid),
      alphaDistributed(processGrid, 1, 1) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityDerivativeEstimationOnOffParallel::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityDerivativeEstimationOnOffParallel::evaluate(DataMatrix& samples,
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

void ModelFittingDensityDerivativeEstimationOnOffParallel::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::fit(DataMatrix& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
  auto& geometryConfig = this->config->getGeometryConfig();
  auto& parallelConfig = this->config->getParallelConfig();
  bool useOfflinePermutation = this->config->getDensityEstimationConfig().useOfflinePermutation_;

  // clear model
  reset();

  // build grid
  gridConfig.dim_ = newDataset.getNcols();
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig, geometryConfig)};

  // build surplus vector
  alpha = DataVector(grid->getSize());

  // Build the offline instance first

  // If the permutation and blow-up approach is applicable to the decomposition type, an object
  // store is given and the offline permutation method is configured, the offline object is obtained
  // from the permutation factory
  if (this->hasObjectStore) {
    const DBMatOffline* objectFromStore =
        this->objectStore->getObject(gridConfig, geometryConfig, refinementConfig,
                                     regularizationConfig, densityEstimationConfig);

    if (objectFromStore != nullptr) {
      offline = std::unique_ptr<DBMatOffline>{objectFromStore->clone()};
    }
  }

  if (DBMatOfflinePermutable::PermutableDecompositions.find(
          densityEstimationConfig.decomposition_) !=
          DBMatOfflinePermutable::PermutableDecompositions.end() &&
      this->hasObjectStore && useOfflinePermutation) {
    // Initialize the permutation factory. If a database path is specified, the path is pased to the
    // permutation factory
    DBMatPermutationFactory permutationFactory;
    if (databaseConfig.filePath_.empty()) {
      permutationFactory = DBMatPermutationFactory(this->objectStore);
    } else {
      permutationFactory = DBMatPermutationFactory(this->objectStore, databaseConfig.filePath_);
    }
    offline = std::unique_ptr<DBMatOffline>{permutationFactory.getPermutedObject(
        config->getGridConfig(), config->getGeometryConfig(), config->getRefinementConfig(),
        config->getRegularizationConfig(), config->getDensityEstimationConfig())};
    offline->interactions = getInteractions(geometryConfig);
  } else if (!databaseConfig.filePath_.empty()) {  // Intialize database if it is provided
    datadriven::DBMatDatabase database(databaseConfig.filePath_);
    // Check if database holds a fitting lhs matrix decomposition
    if (database.hasDataMatrix(gridConfig, refinementConfig, regularizationConfig,
                               densityEstimationConfig)) {
      std::string offlineFilepath = database.getDataMatrix(
          gridConfig, refinementConfig, regularizationConfig, densityEstimationConfig);
      offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildFromFile(offlineFilepath)};
    }
  }

  // Build and decompose offline object if not loaded from database
  if (!offline) {
    // Build offline object by factory, build matrix and decompose
    offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildOfflineObject(
        gridConfig, refinementConfig, regularizationConfig, densityEstimationConfig)};
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
    offline->interactions = getInteractions(geometryConfig);

    if (this->hasObjectStore) {
      this->objectStore->putObject(gridConfig, geometryConfig, refinementConfig,
                                   regularizationConfig, densityEstimationConfig, offline->clone());
    }
  }

  alphaDistributed =
      DataVectorDistributed(processGrid, grid->getSize(), parallelConfig.rowBlockSize_);

  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(
      *offline, *grid, regularizationConfig.lambda_, 0, densityEstimationConfig.decomposition_)};

  online->syncDistributedDecomposition(processGrid, parallelConfig);

  // in SMW decomposition type case, the inverse of the matrix needs to be computed explicitly
  if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
      densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
    offline->compute_inverse_parallel(processGrid, parallelConfig);
  }

  online->computeDensityDerivativeFunctionParallel(
      alphaDistributed, newDataset, *grid, this->config->getDensityEstimationConfig(),
      this->config->getParallelConfig(), processGrid, true,
      this->config->getCrossvalidationConfig().enable_);
  online->setBeta(this->config->getLearnerConfig().forgetRate_);

  alpha = alphaDistributed.toLocalDataVectorBroadcast();

  if (densityEstimationConfig.normalize_) {
    online->normalize(alpha, *grid);
  }
}

bool ModelFittingDensityDerivativeEstimationOnOffParallel::adapt(
    size_t newNoPoints, std::vector<size_t>& deletedGridPoints) {
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

  // Update online object: lhs, rhs and recompute the density function based on the b stored
  online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(), *grid,
                                          newNoPoints - oldNoPoints, deletedGridPoints,
                                          config->getRegularizationConfig().lambda_);
  online->updateRhs(newNoPoints, deletedGridPoints);

  online->syncDistributedDecomposition(processGrid, parallelConfig);
  return true;
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Update the fit (streaming)

    online->computeDensityDerivativeFunctionParallel(
        alphaDistributed, newDataset, *grid, this->config->getDensityEstimationConfig(),
        this->config->getParallelConfig(), processGrid, true,
        this->config->getCrossvalidationConfig().enable_);

    alpha = alphaDistributed.toLocalDataVectorBroadcast();

    if (this->config->getDensityEstimationConfig().normalize_) {
      online->normalize(alpha, *grid);
    }
  }
}

double ModelFittingDensityDerivativeEstimationOnOffParallel::computeResidual(
    DataMatrix& validationData) const {
  DataVectorDistributed bValidation = online->computeWeightedDerivativeBFromBatchParallel(
      validationData, *grid, this->config->getDensityEstimationConfig(),
      this->config->getParallelConfig(), this->processGrid, true);

  DataMatrixDistributed rMatrix = online->getOfflineObject().getUnmodifiedRDistributed(
      this->processGrid, this->config->getParallelConfig());

#ifdef USE_SCALAPACK
  // R * alpha - b_val
  rMatrix.mult(alphaDistributed, bValidation, false, 1.0, -1.0);

  DataVector result = bValidation.toLocalDataVectorBroadcast();
#else
  throw base::not_implemented_exception("built without ScaLAPACK");
#endif /* USE_SCALAPACK */

  return result.l2Norm();
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::updateRegularization(double lambda) {
  if (grid != nullptr) {
    auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
    auto& parallelConfig = this->config->getParallelConfig();

    this->online->getOfflineObject().updateRegularizationParallel(lambda, this->processGrid,
                                                                  parallelConfig);

    // in SMW decomposition type case, the inverse of the matrix needs to be
    // computed explicitly
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
        densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
      online->getOfflineObject().compute_inverse_parallel(processGrid, parallelConfig);
    }
  }
}

bool ModelFittingDensityDerivativeEstimationOnOffParallel::isRefinable() {
  if (grid != nullptr) {
    return online->getOfflineObject().isRefineable();
  }
  return false;
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::reset() {
  grid.reset();
  online.reset();
  refinementsPerformed = 0;
}

void ModelFittingDensityDerivativeEstimationOnOffParallel::resetTraining() {
  if (grid != nullptr) {
    alpha = DataVector(grid->getSize());
    this->online->resetTraining();
  }
}

std::shared_ptr<BlacsProcessGrid>
ModelFittingDensityDerivativeEstimationOnOffParallel::getProcessGrid() const {
  return processGrid;
}

}  // namespace datadriven
}  // namespace sgpp
