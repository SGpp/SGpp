// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>

#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/fista/ElasticNetFunction.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/GroupLassoFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>

#include <cassert>
#include <limits>
#include <random>
#include <set>
#include <vector>

namespace sgpp {
namespace datadriven {

RegressionLearner::RegressionLearner(base::RegularGridConfiguration gridConfig,
                                     base::AdaptivityConfiguration adaptivityConfig,
                                     solver::SLESolverConfiguration solverConfig,
                                     solver::SLESolverConfiguration finalSolverConfig,
                                     datadriven::RegularizationConfiguration regularizationConfig,
                                     std::set<std::set<size_t>> terms)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      finalSolverConfig(finalSolverConfig),
      regularizationConfig(regularizationConfig),
      terms(terms) {
  initializeGrid(gridConfig);
}

RegressionLearner::RegressionLearner(base::RegularGridConfiguration gridConfig,
                                     base::AdaptivityConfiguration adaptivityConfig,
                                     solver::SLESolverConfiguration solverConfig,
                                     solver::SLESolverConfiguration finalSolverConfig,
                                     datadriven::RegularizationConfiguration regularizationConfig)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      finalSolverConfig(finalSolverConfig),
      regularizationConfig(regularizationConfig),
      terms() {
  initializeGrid(gridConfig);
}

void RegressionLearner::train(base::DataMatrix& trainDataset, base::DataVector& classes) {
  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "RegressionLearner::train: length of classes vector does not match to "
        "dataset!");
  }
  auto solver = createSolver(classes.getSize());

  if (solver.type == Solver::solverCategory::cg) {
    systemMatrix = createDMSystem(trainDataset);
  }

  std::unique_ptr<base::OperationMultipleEval> op(
      op_factory::createOperationMultipleEval(*grid, trainDataset));

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
  std::unique_ptr<base::OperationMultipleEval> multOp(
      op_factory::createOperationMultipleEval(*grid, data));
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
      const double L = solver.getL();
      solver.solveFista(*op, weights, classes, solverConfig.maxIterations_, solverConfig.threshold_,
                        L);
      break;
    }
    case Solver::solverCategory::none:
      throw base::application_exception("RegressionLearner::fit: Solver not supported!");
  }
}

void RegressionLearner::refine(base::DataMatrix& data, base::DataVector& classes) {
  // First calculate the training errors for the dataset.
  auto error = predict(data);
  error.sub(classes);
  error.sqr();

  // Calculate the weighted errors per basis function.
  std::unique_ptr<base::OperationMultipleEval> multOp(
      op_factory::createOperationMultipleEval(*grid, data));
  auto errors = base::DataVector(weights.getSize());
  multOp->multTranspose(error, errors);
  errors.componentwise_mult(weights);

  // Refine the grid using the weighted errors.
  auto refineFunctor = base::SurplusRefinementFunctor(errors, adaptivityConfig.numRefinementPoints_,
                                                      adaptivityConfig.refinementThreshold_);
  if (terms.size() > 0) {
    grid->getGenerator().refineInter(refineFunctor, terms);
  } else {
    grid->getGenerator().refine(refineFunctor);
  }

  // tell the SLE manager that the grid changed (for internal
  // data structures)
  if (systemMatrix != nullptr) {
    systemMatrix->prepareGrid();
  }
  weights.resizeZero(grid->getSize());
}

base::Grid& RegressionLearner::getGrid() { return *grid; }

size_t RegressionLearner::getGridSize() const { return grid->getSize(); }

base::DataVector RegressionLearner::getWeights() const { return weights; }

void RegressionLearner::setWeights(base::DataVector weights) { this->weights = weights; }

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

  if (terms.size() > 0) {
    grid->getGenerator().regularInter(gridConfig.level_, terms, gridConfig.t_);
  } else {
    grid->getGenerator().regular(gridConfig.level_, gridConfig.t_);
  }
  weights = base::DataVector(grid->getSize());
}

// maybe pass regularizationConfig instead of state.
std::unique_ptr<datadriven::DMSystemMatrixBase> RegressionLearner::createDMSystem(
    base::DataMatrix& trainDataset) {
  using datadriven::RegularizationType;
  // initializing this is sadly neccesary to resolve a face-off between gcc and clang warnings
  base::OperationMatrix* opMatrix = sgpp::op_factory::createOperationIdentity(*grid);
  switch (regularizationConfig.type_) {
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
    case RegularizationType::GroupLasso:
    case RegularizationType::Lasso:
    case RegularizationType::ElasticNet:
      throw base::application_exception(
          "RegressionLearner::createDMSystem: An unsupported regularization type was chosen!");
  }
  return std::make_unique<datadriven::DMSystemMatrix>(
      *grid, trainDataset, std::shared_ptr<base::OperationMatrix>(opMatrix),
      regularizationConfig.lambda_);
}

RegressionLearner::Solver RegressionLearner::createSolver(size_t n_rows) {
  using solver::SLESolverType;
  switch (solverConfig.type_) {
    case SLESolverType::CG:
      return Solver(std::make_unique<solver::ConjugateGradients>(solverConfig.maxIterations_,
                                                                 solverConfig.eps_));
    case SLESolverType::BiCGSTAB:
      return Solver(
          std::make_unique<solver::BiCGStab>(solverConfig.maxIterations_, solverConfig.eps_));
    case SLESolverType::FISTA:
      return createSolverFista(n_rows);
  }

  throw base::application_exception(
      "RegressionLearner::createSolver: An unsupported solver type was chosen!");
}

RegressionLearner::Solver RegressionLearner::createSolverFista(size_t n_rows) {
  // The FISTA-solver solves loss + lambda * regularization_penalty.
  // We adjust it to align to function like the CG solver, by solving
  // loss + n * lambda * regularization_penalty instead.
  const double lambda = static_cast<double>(n_rows) * regularizationConfig.lambda_;
  using datadriven::RegularizationType;
  switch (regularizationConfig.type_) {
    case RegularizationType::Identity:
      return Solver(
          std::make_unique<solver::Fista<solver::RidgeFunction>>(solver::RidgeFunction(lambda)));
    case RegularizationType::Lasso:
      return Solver(
          std::make_unique<solver::Fista<solver::LassoFunction>>(solver::LassoFunction(lambda)));
    case RegularizationType::ElasticNet:
      return Solver(std::make_unique<solver::Fista<solver::ElasticNetFunction>>(
          solver::ElasticNetFunction(lambda, regularizationConfig.l1Ratio_)));
    case RegularizationType::GroupLasso:
      return Solver(std::make_unique<solver::Fista<solver::GroupLassoFunction>>(
          solver::GroupLassoFunction(lambda, &grid->getStorage())));
    // The following methods are not supported by FISTA.
    case RegularizationType::Diagonal:
    case RegularizationType::Laplace:
      break;
  }

  throw base::application_exception(
      "RegressionLearner::createSolverFista: Regularization type not supported by FISTA!");
}

}  // namespace datadriven
}  // namespace sgpp
