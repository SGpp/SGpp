// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/globaldef.hpp>

// TODO(lettrich): use the system matrix with flexible regularization

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;

using sgpp::base::application_exception;

using sgpp::solver::SLESolver;
using sgpp::solver::ConjugateGradients;
using sgpp::solver::BiCGStab;
using sgpp::solver::SLESolverType;

namespace sgpp {
namespace datadriven {

ModelFittingLeastSquares::ModelFittingLeastSquares(
    std::shared_ptr<DataMiningConfigurationLeastSquares> config)
    : datadriven::ModelFittingBase(), config(config) {}

ModelFittingLeastSquares::~ModelFittingLeastSquares() {}

void ModelFittingLeastSquares::fit(Dataset& dataset) {
  // clear model
  alpha.reset();
  grid.reset();
  systemMatrix.reset();

  // rebuild grid
  initializeGrid(config->getGridConfig());
  // rebuild surplus vector
  alpha = std::make_shared<DataVector>(grid->getSize());

  // create sytem matrix
  systemMatrix = std::shared_ptr<DMSystemMatrixBase>(
      buildSystemMatrix(dataset.getData(), config->getLambda()));

  // create right hand side and system matrix
  auto b = std::make_unique<DataVector>(grid->getSize());
  systemMatrix->generateb(dataset.getTargets(), *b);

  if (config->getSolverRefineConfig().type_ == SLESolverType::CG) {
    solver = std::make_unique<ConjugateGradients>(config->getSolverRefineConfig().maxIterations_,
                                                  config->getSolverRefineConfig().eps_);
  } else if (config->getSolverRefineConfig().type_ == SLESolverType::BiCGSTAB) {
    solver = std::make_unique<BiCGStab>(config->getSolverRefineConfig().maxIterations_,
                                        config->getSolverRefineConfig().eps_);
  } else {
    throw application_exception(
        "LearnerBase::train: An unsupported SLE solver type was "
        "chosen!");
  }

  if (config->getRefinementConfig().numRefinements_ == 0) {
    solver->setMaxIterations(config->getSolverFinalConfig().maxIterations_);
    solver->setEpsilon(config->getSolverFinalConfig().eps_);
  }

  solver->solve(*systemMatrix, *alpha, *b, true, true, DEFAULT_RES_THRESHOLD);
}

void ModelFittingLeastSquares::refine() {
  // create refinement functor
  SurplusRefinementFunctor refinementFunctor(*alpha, config->getRefinementConfig().noPoints_,
                                             config->getRefinementConfig().threshold_);
  // refine grid
  grid->getGenerator().refine(refinementFunctor);

  // tell the SLE manager that the grid changed (for interal data structures)
  // systemMatrix->prepareGrid(); -> empty statement!
  alpha->resizeZero(grid->getSize());

  // if (i == configuration.adaptivityConfig.numRefinements_) {
  solver->setMaxIterations(config->getSolverFinalConfig().maxIterations_);
  solver->setEpsilon(config->getSolverFinalConfig().eps_);
  //}
}

void ModelFittingLeastSquares::update(Dataset& dataset) {
  // create sytem matrix
  systemMatrix.reset();
  systemMatrix = std::shared_ptr<DMSystemMatrixBase>(
      buildSystemMatrix(dataset.getData(), config->getLambda()));

  auto b = std::make_unique<DataVector>(grid->getSize());
  systemMatrix->generateb(dataset.getTargets(), *b);

  solver->solve(*systemMatrix, *alpha, *b, true, true, DEFAULT_RES_THRESHOLD);
}

std::unique_ptr<DMSystemMatrixBase> ModelFittingLeastSquares::buildSystemMatrix(
    DataMatrix& trainDataset, double lambda) {
  auto tmp = std::make_unique<SystemMatrixLeastSquaresIdentity>(*grid, trainDataset, lambda);
  tmp->setImplementation(implementationConfig);
  std::unique_ptr<DMSystemMatrixBase> systemMatrix = std::move(tmp);

  return systemMatrix;
}
}  // namespace datadriven
}  // namespace sgpp
