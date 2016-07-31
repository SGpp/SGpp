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
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/ElasticNetFunction.hpp>
#include <limits>
#include <random>
#include <cassert>

namespace sgpp {
namespace datadriven {

RegressionLearner::RegressionLearner(base::RegularGridConfiguration gridConfig,
                                     base::AdpativityConfiguration adaptivityConfig,
                                     solver::SLESolverConfiguration solverConfig,
                                     solver::SLESolverConfiguration finalSolverConfig,
                                     datadriven::RegularizationConfiguration regularizationConfig)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      finalSolverConfig(finalSolverConfig),
      regularizationConfig(regularizationConfig) {
  initializeGrid(gridConfig);
}

void RegressionLearner::train(base::DataMatrix& trainDataset, base::DataVector& classes) {
  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "RegressionLearner::train: length of classes vector does not match to "
        "dataset!");
  }
  auto solver = std::move(createSolver());

  if (solver.type == Solver::solverCategory::cg) {
    systemMatrix = createDMSystem(trainDataset);
  }

  op = sgpp::op_factory::createOperationMultipleEval(*grid, trainDataset);

  for (size_t curStep = 0; curStep <= adaptivityConfig.numRefinements_; ++curStep) {
    if (curStep > 0) {
      refine(trainDataset, classes);
    }
    if (curStep == adaptivityConfig.numRefinements_) {
      solverConfig = finalSolverConfig;
    }
    fit(solver, classes);
  }
}

base::DataVector RegressionLearner::predict(base::DataMatrix& data) {
  auto prediction = base::DataVector(data.getNrows());
  auto multOp = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  multOp->mult(weights, prediction);
  return prediction;
}

double RegressionLearner::getMSE(base::DataMatrix& data, const base::DataVector& y) {
  const auto yPrediction = predict(data);
  return getMSE(y, yPrediction);
}

void RegressionLearner::fit(Solver& solver, base::DataVector& classes) {
  switch (solver.type) {
    case Solver::solverCategory::cg: {
      auto b = base::DataVector(weights.getSize());
      systemMatrix->generateb(classes, b);
      solver.solveCG(*systemMatrix, weights, b, true, false, solverConfig.threshold_);
      break;
    }
    case Solver::solverCategory::fista: {
      solver.solveFista(*op, weights, classes, solverConfig.maxIterations_,
                        solverConfig.threshold_);
      break;
    }
    case Solver::solverCategory::none:
    default:
      throw base::application_exception("RegressionLearner::fit: Solver not supported!");
  }
}

void RegressionLearner::refine(base::DataMatrix& data, base::DataVector& classes) {
  // First calculate the training errors for the dataset.
  auto error = predict(data);
  error.sub(classes);
  error.sqr();

  // Calculate the weighted errors per basis function.
  auto multOp = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  auto errors = base::DataVector(weights.getSize());
  multOp->multTranspose(error, errors);
  errors.componentwise_mult(weights);

  // Refine the grid using the weighted errors.
  auto refineFunctor = base::SurplusRefinementFunctor(errors, adaptivityConfig.noPoints_,
                                                      adaptivityConfig.threshold_);
  grid->getGenerator().refine(refineFunctor);

  // tell the SLE manager that the grid changed (for internal
  // data structures)
  if (systemMatrix != nullptr) {
    systemMatrix->prepareGrid();
  }
  weights.resizeZero(grid->getSize());
}

size_t RegressionLearner::getGridSize() const { return grid->getSize(); }

void RegressionLearner::initializeWeights() {
  const size_t size = grid->getSize();
  const size_t dimensions = grid->getDimension();

  weights = base::DataVector(size);
  auto& gridStorage = grid->getStorage();

  std::mt19937_64 gen(42);

  const double exponentBase = 0.25;
  for (size_t i = 0; i < size; ++i) {
    base::GridStorage::index_pointer gridIndex = gridStorage.get(i);
    base::GridIndex::level_type levelSum = gridIndex->getLevelSum();
    const double exponent = (static_cast<double>(levelSum) - static_cast<double>(dimensions));
    const double multiplicator = std::pow(exponentBase, exponent);

    std::normal_distribution<> dist(0, std::sqrt(multiplicator));
    weights[i] = dist(gen);
  }
}

base::DataVector RegressionLearner::getWeights() const { return weights; }

double RegressionLearner::getMSE(const base::DataVector& y, base::DataVector yPrediction) {
  yPrediction.sub(y);
  yPrediction.sqr();
  const double totalError = yPrediction.sum();
  return (totalError / static_cast<double>(yPrediction.getSize()));
}

void RegressionLearner::initializeGrid(const base::RegularGridConfiguration gridConfig) {
  using base::GridType;
  // no switch used to avoid missing case warnings
  if (gridConfig.type_ == GridType::LinearBoundary) {
    grid = std::make_unique<base::LinearBoundaryGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::ModLinear) {
    grid = std::make_unique<base::ModLinearGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == GridType::Linear) {
    grid = std::make_unique<base::LinearGrid>(gridConfig.dim_);
  } else {
    throw base::application_exception(
        "RegressionLearner::InitializeGrid: An unsupported grid type was chosen!");
  }

  grid->getGenerator().regular(gridConfig.level_, gridConfig.t_);
  weights = base::DataVector(grid->getSize());
  weights.setAll(0.0);
}

// maybe pass regularizationConfig instead of state.
std::unique_ptr<datadriven::DMSystemMatrixBase> RegressionLearner::createDMSystem(
    base::DataMatrix& trainDataset) {
  using datadriven::RegularizationType;
  std::unique_ptr<base::OperationMatrix> opMatrix;
  switch (regularizationConfig.regType_) {
    case RegularizationType::Identity:
      opMatrix = sgpp::op_factory::createOperationIdentity(*grid);
      break;
    case RegularizationType::Laplace:
      opMatrix = sgpp::op_factory::createOperationLaplace(*grid);
      break;
    case RegularizationType::Diagonal:
      opMatrix =
          sgpp::op_factory::createOperationDiagonal(*grid, regularizationConfig.exponentBase_);
      break;
    case RegularizationType::Lasso:
    case RegularizationType::ElasticNet:
    default:
      throw base::application_exception(
          "RegressionLearner::createDMSystem: An unsupported regularization type was chosen!");
  }
  return std::make_unique<datadriven::DMSystemMatrix>(*grid, trainDataset, std::move(opMatrix),
                                                      regularizationConfig.lambda_);
}

RegressionLearner::Solver RegressionLearner::createSolver() {
  using solver::SLESolverType;
  switch (solverConfig.type_) {
    case SLESolverType::CG:
      return Solver(std::move(std::make_unique<solver::ConjugateGradients>(
          solverConfig.maxIterations_, solverConfig.eps_)));
    case SLESolverType::BiCGSTAB:
      return Solver(std::move(
          std::make_unique<solver::BiCGStab>(solverConfig.maxIterations_, solverConfig.eps_)));
    case SLESolverType::FISTA:
      return createSolverFista();
    default:
      throw base::application_exception(
          "RegressionLearner::createSolver: An unsupported solver type was chosen!");
  }
}

RegressionLearner::Solver RegressionLearner::createSolverFista() {
  using datadriven::RegularizationType;
  switch (regularizationConfig.regType_) {
    case RegularizationType::Identity:
      return Solver(std::make_unique<solver::Fista<solver::RidgeFunction>>(
          solver::RidgeFunction(regularizationConfig.lambda_)));
    case RegularizationType::Lasso:
      return Solver(std::make_unique<solver::Fista<solver::LassoFunction>>(
          solver::LassoFunction(regularizationConfig.lambda_)));
    case RegularizationType::ElasticNet:
      return Solver(std::make_unique<solver::Fista<solver::ElasticNetFunction>>(
          solver::ElasticNetFunction(regularizationConfig.lambda_, regularizationConfig.l1Ratio_)));
    // The following methods are not supported by FISTA.
    case RegularizationType::Diagonal:
    case RegularizationType::Laplace:
    default:
      throw base::application_exception(
          "RegressionLearner::createSolverFista: Regularization type not supported by FISTA!");
  }
}

}  // namespace datadriven
}  // namespace sgpp
