// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

// TODO(lettrich): use the system matrix with flexible regularization

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

ModelFittingLeastSquares::ModelFittingLeastSquares(
    sgpp::datadriven::DataMiningConfigurationLeastSquares config)
    : datadriven::ModelFittingBase(), configuration(config) {}

ModelFittingLeastSquares::~ModelFittingLeastSquares() {}

datadriven::DMSystemMatrixBase* ModelFittingLeastSquares::createSystemMatrix(
    base::DataMatrix& trainDataset, double lambda) {
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

  initializeGrid(configuration.getGridConfig());
  alpha = std::make_shared<DataVector>(grid->getSize());
  alpha->setAll(0);

  // create DMSystem
  systemMatrix = std::shared_ptr<datadriven::DMSystemMatrixBase>(
      createSystemMatrix(dataset.getData(), configuration.getLambda()));

  if (configuration.getSolverRefineConfig().type_ == sgpp::solver::SLESolverType::CG) {
    solver = std::make_shared<solver::ConjugateGradients>(
        configuration.getSolverRefineConfig().maxIterations_,
        configuration.getSolverRefineConfig().eps_);
  } else if (configuration.getSolverRefineConfig().type_ == sgpp::solver::SLESolverType::BiCGSTAB) {
    solver =
        std::make_shared<solver::BiCGStab>(configuration.getSolverRefineConfig().maxIterations_,
                                           configuration.getSolverRefineConfig().eps_);
  } else {
    throw base::application_exception(
        "LearnerBase::train: An unsupported SLE solver type was "
        "chosen!");
  }

  DataVector b(grid->getSize());
  systemMatrix->generateb(dataset.getTargets(), b);

  if (configuration.getRefinementConfig().numRefinements_ == 0) {
    solver->setMaxIterations(configuration.getSolverFinalConfig().maxIterations_);
    solver->setEpsilon(configuration.getSolverFinalConfig().eps_);
  }

  solver->solve(*systemMatrix, *alpha, b, true, true, DEFAULT_RES_THRESHOLD);

  //  double tmp1, tmp2, tmp3, tmp4;
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
}  // namespace datadriven
}  // namespace sgpp
