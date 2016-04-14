// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/RegressionLearner.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>

namespace sgpp {
namespace datadriven {

RegressionLearner::RegressionLearner(
    sgpp::base::RegularGridConfiguration gridConfig,
    sgpp::base::AdpativityConfiguration adaptivityConfig,
    sgpp::solver::SLESolverConfiguration solverConfig,
    sgpp::datadriven::RegularizationConfiguration regularizationConfig)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig) {
  // initialize grid, etc.
  initializeGrid(gridConfig);  // do I need to store this?
}

void RegressionLearner::train(sgpp::base::DataMatrix& trainDataset,
                              sgpp::base::DataVector& classes) {
  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "RegressionLearner::train: length of classes vector does not match to "
        "dataset!");
  }
  const auto lambda = 0.001;  // TODO(krenzls): pass regularization param.
  const auto DMSystem = createDMSystem(trainDataset, lambda);
  std::unique_ptr<sgpp::solver::SLESolver> solver;
  // TODO(krenzls): support multiple solver types!
  solver = std::make_unique<sgpp::solver::ConjugateGradients>(solverConfig.maxIterations_,
                                                              solverConfig.eps_);

  for (size_t i = 0; i < adaptivityConfig.numRefinements_ + 1; ++i) {
    // We don't need to refine the grid during the first iteration!
    if (i > 0) {
      auto refineFunctor = sgpp::base::SurplusRefinementFunctor(weights, adaptivityConfig.noPoints_,
                                                                adaptivityConfig.threshold_);
      grid->getGenerator().refine(refineFunctor);

      // tell the SLE manager that the grid changed (for internal
      // data structures)
      DMSystem->prepareGrid();
      weights.resizeZero(grid->getSize());
    } else {
      // no refinement needed!
    }

    auto b = sgpp::base::DataVector(weights.getSize());
    DMSystem->generateb(classes, b);

    if (i == adaptivityConfig.numRefinements_) {
      // change the configuration of the solver for this last step
      solver->setMaxIterations(solverConfig.maxIterations_);
      solver->setEpsilon(solverConfig.eps_);
    }

    solver->solve(*DMSystem, weights, b, true, false, 0.0);
  }
}

void RegressionLearner::initializeGrid(const sgpp::base::RegularGridConfiguration gridConfig) {
  using sgpp::base::GridType;
  // no switch used to avoid missing case warnings
  if (gridConfig.type_ == GridType::LinearBoundary) {
    grid = std::make_unique<sgpp::base::LinearBoundaryGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::ModLinear) {
    grid = std::make_unique<sgpp::base::ModLinearGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::Linear) {
    grid = std::make_unique<sgpp::base::LinearGrid>(gridConfig.dim_);
  } else {
    throw base::application_exception(
        "RegressionLearner::InitializeGrid: An unsupported grid type was chosen!");
  }

  grid->getGenerator().regular(gridConfig.level_);
  weights = sgpp::base::DataVector(grid->getSize());
  weights.setAll(0.0);
}

// maybe pass regularizationConfig instead of state.
std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> RegressionLearner::createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda) {
  using datadriven::RegularizationType;
  switch (regularizationConfig.regType_) {
    case RegularizationType::Identity:
      opMatrix = sgpp::op_factory::createOperationIdentity(*grid);
      break;
    case RegularizationType::Laplace:
      opMatrix = sgpp::op_factory::createOperationLaplace(*grid);
      break;
    default:
      throw base::application_exception(
          "RegressionLearner::createDMSystem: An unsupported regularization type was chosen!");
  }
  return std::make_unique<sgpp::datadriven::DMSystemMatrix>(*grid, trainDataset, *opMatrix, lambda);
}

}  // namespace datadriven
}  // namespace sgpp
