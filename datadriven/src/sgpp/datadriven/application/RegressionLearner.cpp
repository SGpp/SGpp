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

void RegressionLearner::train() {
  // do something ;)
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
  std::unique_ptr<sgpp::base::OperationMatrix> opMatrix;
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
