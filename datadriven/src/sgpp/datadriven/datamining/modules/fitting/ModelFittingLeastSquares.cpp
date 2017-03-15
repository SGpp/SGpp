// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
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

ModelFittingLeastSquares::ModelFittingLeastSquares(const FitterConfigurationLeastSquares& config)
    : ModelFittingBase{}, systemMatrix{nullptr} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationLeastSquares>(config));
  solver = std::unique_ptr<SLESolver>{buildSolver(this->config->getSolverFinalConfig())};
}

double ModelFittingLeastSquares::evaluate(const DataVector& sample) const {
  auto opEval = std::unique_ptr<base::OperationEval>{op_factory::createOperationEval(*grid)};
  return opEval->eval(alpha, sample);
}

void ModelFittingLeastSquares::evaluate(DataMatrix& samples, DataVector& results) {
  auto opMultEval = std::unique_ptr<base::OperationMultipleEval>{
      op_factory::createOperationMultipleEval(*grid, samples, config->getMultipleEvalConfig())};
  opMultEval->eval(alpha, results);
}

void ModelFittingLeastSquares::fit(Dataset& dataset) {
  // clear model
  grid.reset();
  systemMatrix.reset();

  // build grid
  auto& gridConfig = config->getGridConfig();
  gridConfig.dim_ = dataset.getDimension();
  grid = std::unique_ptr<Grid>{buildGrid(config->getGridConfig())};
  // build surplus vector
  alpha = DataVector(grid->getSize());

  // create sytem matrix
  systemMatrix = std::unique_ptr<DMSystemMatrixBase>(
      buildSystemMatrix(*grid, dataset.getData(), config->getRegularizationConfig().lambda_,
                        config->getMultipleEvalConfig()));

  // create right hand side and system matrix
  auto b = std::make_unique<DataVector>(grid->getSize());
  systemMatrix->generateb(dataset.getTargets(), *b);

  reconfigureSolver(*solver, config->getSolverFinalConfig());
  solver->solve(*systemMatrix, alpha, *b, true, true, DEFAULT_RES_THRESHOLD);
}

void ModelFittingLeastSquares::refine() {
  if (grid != nullptr) {
    // create refinement functor
    SurplusRefinementFunctor refinementFunctor(alpha, config->getRefinementConfig().noPoints_,
                                               config->getRefinementConfig().threshold_);
    // refine grid
    grid->getGenerator().refine(refinementFunctor);

    // tell the SLE manager that the grid changed (for interal data structures)
    // systemMatrix->prepareGrid(); -> empty statement!
    alpha.resizeZero(grid->getSize());

    throw application_exception(
        "ModelFittingLeastSquares: Can't refine before initial grid is created");
  }
}

void ModelFittingLeastSquares::update(Dataset& dataset) {
  if (grid != nullptr) {
    // create sytem matrix
    systemMatrix.reset();
    systemMatrix = std::unique_ptr<DMSystemMatrixBase>(
        buildSystemMatrix(*grid, dataset.getData(), config->getRegularizationConfig().lambda_,
                          config->getMultipleEvalConfig()));

    auto b = std::make_unique<DataVector>(grid->getSize());
    systemMatrix->generateb(dataset.getTargets(), *b);

    reconfigureSolver(*solver, config->getSolverRefineConfig());
    solver->solve(*systemMatrix, alpha, *b, true, true, DEFAULT_RES_THRESHOLD);
  } else {
    fit(dataset);
  }
}

DMSystemMatrixBase* ModelFittingLeastSquares::buildSystemMatrix(
    Grid& grid, DataMatrix& trainDataset, double lambda,
    OperationMultipleEvalConfiguration& mutipleEvalconfig) const {
  auto systemMatrix = new SystemMatrixLeastSquaresIdentity(grid, trainDataset, lambda);
  systemMatrix->setImplementation(mutipleEvalconfig);

  return systemMatrix;
}

}  // namespace datadriven
}  // namespace sgpp
