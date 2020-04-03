// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixDensityRatioEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityRatioEstimation.hpp>
#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

// TODO(lettrich): allow different refinement types
// TODO(lettrich): allow different refinement criteria

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;

using sgpp::base::application_exception;

using sgpp::solver::SLESolver;

namespace sgpp {
namespace datadriven {

ModelFittingDensityRatioEstimation::ModelFittingDensityRatioEstimation(
    const FitterConfigurationLeastSquares &config)
    : ModelFittingBaseSingleGrid{}, refinementsPerformed{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationLeastSquares>(config));
  solver = std::unique_ptr<SLESolver>{buildSolver(this->config->getSolverFinalConfig())};
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityRatioEstimation::evaluate(const DataVector &sample) {
  auto opEval = std::unique_ptr<base::OperationEval>{op_factory::createOperationEval(*grid)};
  return opEval->eval(alpha, sample);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityRatioEstimation::evaluate(DataMatrix &samples, DataVector &results) {
  auto opMultEval = std::unique_ptr<base::OperationMultipleEval>{
      op_factory::createOperationMultipleEval(*grid, samples, config->getMultipleEvalConfig())};
  opMultEval->eval(alpha, results);
}

void ModelFittingDensityRatioEstimation::fit(Dataset &newDatasetP, Dataset &newDatasetQ) {
  // clear model
  reset();
  this->dataset = &newDatasetP;
  this->extraDataset = &newDatasetQ;

  // build grid
  auto &gridConfig = config->getGridConfig();
  gridConfig.dim_ = dataset->getDimension();  // extraDataset works as well
  grid = std::unique_ptr<Grid>{buildGrid(config->getGridConfig())};
  // build surplus vector
  alpha = DataVector(grid->getSize());

  assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
}

bool ModelFittingDensityRatioEstimation::adapt() {
  if (grid != nullptr) {
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      SurplusRefinementFunctor refinementFunctor(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
      // refine grid
      auto noPoints = grid->getSize();

      grid->getGenerator().refine(refinementFunctor);

      if (grid->getSize() > noPoints) {
        // Tell the SLE manager that the grid changed (for internal data structures)
        alpha.resizeZero(grid->getSize());

        assembleSystemAndSolve(config->getSolverRefineConfig(), alpha);
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
        "ModelFittingDensityRatioEstimation: Can't refine before initial grid is created");
  }
}

void ModelFittingDensityRatioEstimation::update(Dataset &newDatasetP, Dataset &newDatasetQ) {
  if (grid != nullptr) {
    reset();
    // reassign datasets
    dataset = &newDatasetP;
    extraDataset = &newDatasetQ;
    // create sytem matrix
    assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
  } else {
    fit(newDatasetP, newDatasetQ);
  }
}

DMSystemMatrixDRE *ModelFittingDensityRatioEstimation::buildSystemMatrix(
    Grid &grid, DataMatrix &trainDatasetP, DataMatrix &trainDatasetQ, double lambda,
    OperationMultipleEvalConfiguration &mutipleEvalconfig) const {
  auto systemMatrix =
      new SystemMatrixDensityRatioEstimation(grid, trainDatasetP, trainDatasetQ, lambda);
  systemMatrix->setImplementation(mutipleEvalconfig);

  return systemMatrix;
}

void ModelFittingDensityRatioEstimation::reset() {
  grid.reset();
  refinementsPerformed = 0;
}

void ModelFittingDensityRatioEstimation::assembleSystemAndSolve(
    const SLESolverConfiguration &solverConfig, DataVector &alpha) const {
  auto systemMatrix = std::unique_ptr<DMSystemMatrixDRE>(buildSystemMatrix(
      *grid, dataset->getData(), extraDataset->getData(), config->getRegularizationConfig().lambda_,
      config->getMultipleEvalConfig()));

  DataVector b(grid->getSize());
  systemMatrix->generateb(b);

  reconfigureSolver(*solver, solverConfig);
  solver->solve(*systemMatrix, alpha, b, true, verboseSolver, DEFAULT_RES_THRESHOLD);
}
}  // namespace datadriven
}  // namespace sgpp
