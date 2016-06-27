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

#include <limits>
#include <random>

namespace sgpp {
namespace datadriven {

RegressionLearner::RegressionLearner(
    sgpp::base::RegularGridConfiguration gridConfig,
    sgpp::base::AdpativityConfiguration adaptivityConfig,
    sgpp::solver::SLESolverConfiguration solverConfig,
    sgpp::solver::SLESolverConfiguration finalSolverConfig,
    sgpp::datadriven::RegularizationConfiguration regularizationConfig)
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      finalSolverConfig(finalSolverConfig),
      regularizationConfig(regularizationConfig) {
  initializeGrid(gridConfig);
}

void RegressionLearner::train(sgpp::base::DataMatrix& trainDataset,
                              sgpp::base::DataVector& classes) {
  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "RegressionLearner::train: length of classes vector does not match to "
        "dataset!");
  }
  const auto DMSystem = createDMSystem(trainDataset);



  for (size_t curStep = 0; curStep <= adaptivityConfig.numRefinements_; ++curStep) {
    if (curStep > 0) {
      refine(*DMSystem, trainDataset, classes);
    }
    if (curStep == adaptivityConfig.numRefinements_) {
      solverConfig = finalSolverConfig; // TODO(krenzls): Add final solver!
    }
    auto solver = createSolver();
    fit(*DMSystem, *solver, classes);
  }
}


sgpp::base::DataVector RegressionLearner::predict(sgpp::base::DataMatrix& data) {
  auto prediction = sgpp::base::DataVector(data.getNrows());
  auto multOp = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  multOp->mult(weights, prediction);
  return prediction;
}

double RegressionLearner::getMSE(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y) {
  const auto yPrediction = predict(data);
  return getMSE(y, yPrediction);
}

void RegressionLearner::fit(sgpp::datadriven::DMSystemMatrixBase& DMSystem,
                            sgpp::solver::SLESolver& solver, sgpp::base::DataVector& classes) {
  auto b = sgpp::base::DataVector(weights.getSize());
  DMSystem.generateb(classes, b);
  solver.solve(DMSystem, weights, b, true, false, 0.0);
}

void RegressionLearner::refine(sgpp::datadriven::DMSystemMatrixBase& DMSystem, sgpp::base::DataMatrix& data,
                               sgpp::base::DataVector& classes) {
  // First calculate the training errors for the dataset.
  auto error = predict(data);
  error.sub(classes);
  error.sqr();

  // Calculate the weighted errors per basis function.
  auto multOp = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  auto errors = sgpp::base::DataVector(weights.getSize());
  multOp->multTranspose(error, errors);
  errors.componentwise_mult(weights);

  // Refine the grid using the weighted errors.
  auto refineFunctor = sgpp::base::SurplusRefinementFunctor(errors, adaptivityConfig.noPoints_,
                                                            adaptivityConfig.threshold_);
  grid->getGenerator().refine(refineFunctor);

  // tell the SLE manager that the grid changed (for internal
  // data structures)
  DMSystem.prepareGrid();
  weights.resizeZero(grid->getSize());
}

size_t RegressionLearner::getGridSize() const {
    return grid->getSize();
}

void RegressionLearner::initializeWeights() {
    const size_t size = grid->getSize();
    const size_t dimensions = grid->getDimension();

    weights = sgpp::base::DataVector(size);
    auto& gridStorage = grid->getStorage();

    std::mt19937_64 gen(42);

    const double exponentBase = 0.25;
    for(size_t i = 0; i < size; ++i) {

        sgpp::base::GridStorage::index_pointer gridIndex = gridStorage.get(i);
        sgpp::base::GridIndex::level_type levelSum = gridIndex->getLevelSum();
        const double exponent = (static_cast<double>(levelSum) - static_cast<double>(dimensions));
        const double multiplicator = std::pow(exponentBase, exponent);

        std::normal_distribution<> dist(0, multiplicator);
        weights[i] = dist(gen);
    }
}

double RegressionLearner::getMSE(const sgpp::base::DataVector& y,
                                 sgpp::base::DataVector yPrediction) {
  yPrediction.sub(y);
  yPrediction.sqr();
  const double totalError = yPrediction.sum();
  return (totalError / static_cast<double>(yPrediction.getSize()));
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

  grid->getGenerator().regular(gridConfig.level_, gridConfig.t_);
  weights = sgpp::base::DataVector(grid->getSize());
  weights.setAll(0.0);
}

// maybe pass regularizationConfig instead of state.
std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> RegressionLearner::createDMSystem(
    sgpp::base::DataMatrix& trainDataset) {
  using datadriven::RegularizationType;
  std::unique_ptr<sgpp::base::OperationMatrix> opMatrix;
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
    default:
      throw base::application_exception(
          "RegressionLearner::createDMSystem: An unsupported regularization type was chosen!");
  }
  return std::make_unique<sgpp::datadriven::DMSystemMatrix>(
      *grid, trainDataset, std::move(opMatrix), regularizationConfig.lambda_);
}

std::unique_ptr<sgpp::solver::SLESolver> RegressionLearner::createSolver() {
  using sgpp::solver::SLESolverType;
  decltype(createSolver()) solver;
  switch (solverConfig.type_) {
    case SLESolverType::CG:
      solver = std::make_unique<sgpp::solver::ConjugateGradients>(solverConfig.maxIterations_,
                                                                  solverConfig.eps_);
      break;
    case SLESolverType::BiCGSTAB:
      solver =
          std::make_unique<sgpp::solver::BiCGStab>(solverConfig.maxIterations_, solverConfig.eps_);
      break;
    default:
      throw base::application_exception(
          "RegressionLearner::createSolver: An unsupported solver type was chosen!");
  }
  return solver;
}

}  // namespace datadriven
}  // namespace sgpp
