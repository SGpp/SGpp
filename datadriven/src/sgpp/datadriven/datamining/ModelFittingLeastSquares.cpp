// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "ModelFittingLeastSquares.hpp"

#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>

// TODO: use the system matrix with flexible regularization

namespace SGPP {
namespace datadriven {

ModelFittingLeastSquares::ModelFittingLeastSquares(
    SGPP::datadriven::DataMiningConfigurationLeastSquares config)
    : datadriven::ModelFittingBase(), configuration(config) {}

ModelFittingLeastSquares::~ModelFittingLeastSquares() {}

datadriven::DMSystemMatrixBase* ModelFittingLeastSquares::createSystemMatrix(
    base::DataMatrix& trainDataset, float_t lambda) {
  datadriven::SystemMatrixLeastSquaresIdentity* systemMatrix =
      new datadriven::SystemMatrixLeastSquaresIdentity(*(this->grid), trainDataset, lambda);
  systemMatrix->setImplementation(this->implementationConfiguration);
  return systemMatrix;
}

void ModelFittingLeastSquares::fit(datadriven::Dataset& dataset) {
  //  LearnerTiming result;
  //  result.timeComplete_ = 0.0;
  //  result.timeMultComplete_ = 0.0;
  //  result.timeMultCompute_ = 0.0;
  //  result.timeMultTransComplete_ = 0.0;
  //  result.timeMultTransCompute_ = 0.0;
  //  result.timeRegularization_ = 0.0;
  //  result.GFlop_ = 0.0;
  //  result.GByte_ = 0.0;

  alpha.reset();
  grid.reset();
  systemMatrix.reset();

  //  InitializeGrid(GridConfig);

  // create DMSystem
  systemMatrix = std::shared_ptr<datadriven::DMSystemMatrixBase>(
      createSystemMatrix(dataset.getData(), configuration.getLambda()));

  if (configuration.getSolverRefineConfig().type_ == SGPP::solver::SLESolverType::CG) {
    solver = std::make_shared<solver::ConjugateGradients>(
        configuration.getSolverRefineConfig().maxIterations_,
        configuration.getSolverRefineConfig().eps_);
  } else if (configuration.getSolverRefineConfig().type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
    solver =
        std::make_shared<solver::BiCGStab>(configuration.getSolverRefineConfig().maxIterations_,
                                           configuration.getSolverRefineConfig().eps_);
  } else {
    throw base::application_exception(
        "LearnerBase::train: An unsupported SLE solver type was "
        "chosen!");
  }

  SGPP::base::DataVector b(alpha->getSize());
  systemMatrix->generateb(dataset.getTargets(), b);

  if (configuration.getRefinementConfig().numRefinements_ == 0) {
    solver->setMaxIterations(configuration.getSolverFinalConfig().maxIterations_);
    solver->setEpsilon(configuration.getSolverFinalConfig().eps_);
  }

  solver->solve(*systemMatrix, *alpha, b, true, false, 0.0);

  //  float_t tmp1, tmp2, tmp3, tmp4;
  //  systemMatrix->getTimers(tmp1, tmp2, tmp3, tmp4);
  //  result.timeComplete_ = execTime_;
  //  result.timeMultComplete_ = tmp1;
  //  result.timeMultCompute_ = tmp2;
  //  result.timeMultTransComplete_ = tmp3;
  //  result.timeMultTransCompute_ = tmp4;
  //  result.timeRegularization_ = 0.0;
  //  result.GFlop_ = GFlop_;
  //  result.GByte_ = GByte_;
}

void ModelFittingLeastSquares::refine() {
  //  // disable refinement here!
  //  auto refinementFunctor =
  //  std::make_shared<base::SurplusRefinementFunctor>(alpha,
  //          configuration.adaptivityConfig.noPoints_,
  //          configuration.adaptivityConfig.threshold_);
  //  std::unique_ptr<base::GridGenerator>(grid->createGridGenerator())->refine(
  //          refinementFunctor.get());
  //
  //  // tell the SLE manager that the grid changed (for interal data
  //  structures)
  //  systemMatrix->prepareGrid();
  //
  //  alpha->resizeZero(grid->getSize());
  //
  //  if (i == configuration.adaptivityConfig.numRefinements_) {
  //    myCG->setMaxIterations(configuration.getSolverFinalConfig().maxIterations_);
  //    myCG->setEpsilon(configuration.getSolverFinalConfig().eps_);
  //  }
}

void ModelFittingLeastSquares::update(datadriven::Dataset& dataset) {
  //  // disable refinement here!
  //  auto refinementFunctor =
  //  std::make_shared<base::SurplusRefinementFunctor>(alpha,
  //          configuration.adaptivityConfig.noPoints_,
  //          configuration.adaptivityConfig.threshold_);
  //  std::unique_ptr<base::GridGenerator>(grid->createGridGenerator())->refine(
  //          refinementFunctor.get());
  //
  //  // tell the SLE manager that the grid changed (for interal data
  //  structures)
  //  systemMatrix->prepareGrid();
  //
  //  alpha->resizeZero(grid->getSize());
  //
  //  if (i == configuration.adaptivityConfig.numRefinements_) {
  //    myCG->setMaxIterations(configuration.getSolverFinalConfig().maxIterations_);
  //    myCG->setEpsilon(configuration.getSolverFinalConfig().eps_);
  //  }
}
}
}
