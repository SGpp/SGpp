// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif /* USE_GSL */

#include <fstream>
#include <iostream>
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

ModelFittingDensityEstimationOnOff::ModelFittingDensityEstimationOnOff(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation() {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  this->hasObjectStore = false;
}

ModelFittingDensityEstimationOnOff::ModelFittingDensityEstimationOnOff(
    const FitterConfigurationDensityEstimation& config,
    std::shared_ptr<DBMatObjectStore> objectStore)
    : ModelFittingDensityEstimationOnOff(config) {
  this->objectStore = objectStore;
  this->hasObjectStore = true;
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityEstimationOnOff::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityEstimationOnOff::evaluate(DataMatrix& samples, DataVector& results) {
  online->eval(alpha, samples, results, *grid);
}

void ModelFittingDensityEstimationOnOff::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationOnOff::fit(DataMatrix& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
  auto& geometryConfig = this->config->getGeometryConfig();
  bool useOfflinePermutation = this->config->getDensityEstimationConfig().useOfflinePermutation_;

  // clear model
  reset();

  // build grid
  gridConfig.dim_ = newDataset.getNcols();
  // std::cout << "Dataset dimension " << gridConfig.dim_ << std::endl;
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
    offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    offline->interactions = getInteractions(geometryConfig);

    if (this->hasObjectStore) {
      this->objectStore->putObject(gridConfig, geometryConfig, refinementConfig,
                                   regularizationConfig, densityEstimationConfig, offline->clone());
    }
  }

  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(
      *offline, *grid, regularizationConfig.lambda_, 0, densityEstimationConfig.decomposition_)};

  // in SMW decomposition type case, the inverse of the matrix needs to be computed explicitly
  if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
      densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
    offline->compute_inverse();
  }

  online->computeDensityFunction(alpha, newDataset, *grid,
                                 this->config->getDensityEstimationConfig(), true,
                                 this->config->getCrossvalidationConfig().enable_);
  online->setBeta(this->config->getLearnerConfig().learningRate_);

  if (densityEstimationConfig.normalize_) {
    online->normalize(alpha, *grid);
  }
}

bool ModelFittingDensityEstimationOnOff::adapt(size_t newNoPoints,
                                               std::vector<size_t>& deletedGridPoints) {
  // Coarsening, remove idx from alpha
  if (deletedGridPoints.size() > 0) {
    // Restructure alpha
    alpha.remove(deletedGridPoints);
  }
  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // Refinement, expand alpha
  if (newNoPoints > oldNoPoints) {
    alpha.resizeZero(newNoPoints);
  }

  // Update online object: lhs, rhs and recompute the density function based on the b stored
  online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(), *grid,
                                          newNoPoints - oldNoPoints, deletedGridPoints,
                                          config->getRegularizationConfig().lambda_);
  online->updateRhs(newNoPoints, deletedGridPoints);
  return true;
}

void ModelFittingDensityEstimationOnOff::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityEstimationOnOff::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Update the fit (streaming)

    online->computeDensityFunction(alpha, newDataset, *grid,
                                   this->config->getDensityEstimationConfig(), true,
                                   this->config->getCrossvalidationConfig().enable_);

    if (this->config->getDensityEstimationConfig().normalize_) {
      online->normalize(alpha, *grid);
    }
  }
}

double ModelFittingDensityEstimationOnOff::computeResidual(DataMatrix& validationData) const {
  DataVector bValidation = online->computeWeightedBFromBatch(
      validationData, *grid, this->config->getDensityEstimationConfig(), true);

  DataMatrix rMatrix = online->getOfflineObject().getUnmodifiedR();

#ifdef USE_GSL

  gsl_matrix_view R_view =
      gsl_matrix_view_array(rMatrix.getPointer(), rMatrix.getNrows(), rMatrix.getNcols());
  gsl_vector_const_view alpha_view =
      gsl_vector_const_view_array(alpha.getPointer(), alpha.getSize());
  gsl_vector_view b_view = gsl_vector_view_array(bValidation.getPointer(), bValidation.getSize());

  // R * alpha - b_val
  gsl_blas_dgemv(CblasNoTrans, 1.0, &R_view.matrix, &alpha_view.vector, -1.0, &b_view.vector);
#else
  throw base::not_implemented_exception("built withot GSL");
#endif /* USE_GSL */

  return bValidation.l2Norm();
}

void ModelFittingDensityEstimationOnOff::updateRegularization(double lambda) {
  if (grid != nullptr) {
    auto& densityEstimationConfig = this->config->getDensityEstimationConfig();

    this->online->getOfflineObject().updateRegularization(lambda);

    // in SMW decomposition type case, the inverse of the matrix needs to be
    // computed explicitly
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
        densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
      online->getOfflineObject().compute_inverse();
    }
  }
}

bool ModelFittingDensityEstimationOnOff::isRefinable() {
  if (grid != nullptr) {
    return online->getOfflineObject().isRefineable();
  }
  return false;
}

void ModelFittingDensityEstimationOnOff::reset() {
  grid.reset();
  online.reset();
  refinementsPerformed = 0;
}

void ModelFittingDensityEstimationOnOff::resetTraining() {
  if (grid != nullptr) {
    alpha = DataVector(grid->getSize());
    this->online->resetTraining();
  }
}

}  // namespace datadriven
}  // namespace sgpp
