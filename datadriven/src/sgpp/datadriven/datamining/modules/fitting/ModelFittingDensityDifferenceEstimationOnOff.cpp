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
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationOnOff.hpp>

#include <string>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::RefinementFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;
using sgpp::base::RefinementFunctorType;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityDifferenceEstimationOnOff::ModelFittingDensityDifferenceEstimationOnOff(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation() {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityDifferenceEstimationOnOff::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityDifferenceEstimationOnOff::evaluate(DataMatrix& samples,
                                                            DataVector& results) {
  online->eval(alpha, samples, results, *grid);
}

void ModelFittingDensityDifferenceEstimationOnOff::fit(Dataset& newDatasetP, Dataset& newDatasetQ) {
  dataset = &newDatasetP;
  extraDataset = &newDatasetQ;
  fit(newDatasetP.getData(), newDatasetQ.getData());
}

void ModelFittingDensityDifferenceEstimationOnOff::fit(DataMatrix& newDatasetP,
                                                       DataMatrix& newDatasetQ) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();

  // clear model
  reset();

  // build grid
  gridConfig.dim_ = newDatasetP.getNcols();  // newDatasetQ.getNcols() works as well
  std::cout << "Dataset dimension " << gridConfig.dim_ << std::endl;
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  // grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  // build surplus vector
  alpha = DataVector(grid->getSize());

  // Intialize database if it is provided
  if (!databaseConfig.filePath_.empty()) {
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
  }

  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(
      *offline, *grid, regularizationConfig.lambda_, 0, densityEstimationConfig.decomposition_)};

  // in SMW decomposition type case, the inverse of the matrix needs to be computed explicitly
  if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_ortho ||
      densityEstimationConfig.decomposition_ == MatrixDecompositionType::SMW_chol) {
    offline->compute_inverse();
  }

  online->computeDensityDifferenceFunction(alpha, newDatasetP, newDatasetQ, *grid,
                                           this->config->getDensityEstimationConfig(), true,
                                           this->config->getCrossvalidationConfig().enable_);
  online->setBeta(this->config->getLearnerConfig().learningRate_);

  if (densityEstimationConfig.normalize_) {
    online->normalize(alpha, *grid);
  }
}

bool ModelFittingDensityDifferenceEstimationOnOff::adapt(size_t newNoPoints,
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

void ModelFittingDensityDifferenceEstimationOnOff::update(Dataset& newDatasetP,
                                                          Dataset& newDatasetQ) {
  dataset = &newDatasetP;
  extraDataset = &newDatasetQ;
  update(newDatasetP.getData(), newDatasetQ.getData());
}

void ModelFittingDensityDifferenceEstimationOnOff::update(DataMatrix& newDatasetP,
                                                          DataMatrix& newDatasetQ) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDatasetP, newDatasetQ);
  } else {
    // Update the fit (streaming)
    online->computeDensityDifferenceFunction(alpha, newDatasetP, newDatasetQ, *grid,
                                             this->config->getDensityEstimationConfig(), true,
                                             this->config->getCrossvalidationConfig().enable_);

    if (this->config->getDensityEstimationConfig().normalize_) {
      online->normalize(alpha, *grid);
    }
  }
}

bool ModelFittingDensityDifferenceEstimationOnOff::isRefinable() {
  if (grid != nullptr) {
    return online->getOfflineObject().isRefineable();
  }
  return false;
}

void ModelFittingDensityDifferenceEstimationOnOff::reset() {
  grid.reset();
  online.reset();
  refinementsPerformed = 0;
}

}  // namespace datadriven
}  // namespace sgpp
