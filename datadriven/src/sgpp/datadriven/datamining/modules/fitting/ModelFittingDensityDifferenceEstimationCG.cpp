// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DensityDifferenceSystemMatrix.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationCG.hpp>

#include <string>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::RefinementFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;
using sgpp::base::RefinementFunctorType;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityDifferenceEstimationCG::ModelFittingDensityDifferenceEstimationCG(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation{}, bNumP{0}, bDenomP{0}, bNumQ{0}, bDenomQ{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityDifferenceEstimationCG::evaluate(const DataVector& sample) {
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
  return opEval->eval(alpha, sample);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityDifferenceEstimationCG::evaluate(DataMatrix& samples, DataVector& results) {
  sgpp::op_factory::createOperationMultipleEval(*grid, samples)->eval(alpha, results);
}

void ModelFittingDensityDifferenceEstimationCG::fit(Dataset& newDatasetP, Dataset& newDatasetQ) {
  dataset = &newDatasetP;
  extraDataset = &newDatasetQ;
  fit(newDatasetP.getData(), newDatasetQ.getData());
}

void ModelFittingDensityDifferenceEstimationCG::fit(DataMatrix& newDatasetP,
                                                    DataMatrix& newDatasetQ) {
  // Reset grid
  reset();

  // Setup new grid
  auto& gridConfig = this->config->getGridConfig();
  gridConfig.dim_ = newDatasetP.getNcols();  // newDatasetQ.getNcols() works as well
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  // build surplus vector
  alpha = DataVector(grid->getSize());

  // Initialize the right hand side (numerator and denominator)
  bNumP = DataVector(grid->getSize());
  bDenomP = DataVector(grid->getSize());
  bNumQ = DataVector(grid->getSize());
  bDenomQ = DataVector(grid->getSize());

  // Now that everythin is set up with zero data, simply call the update method
  update(newDatasetP, newDatasetQ);
}

bool ModelFittingDensityDifferenceEstimationCG::adapt(size_t newNoPoints,
                                                      std::vector<size_t>& deletedGridPoints) {
  // Coarsening, remove idx from alpha
  if (deletedGridPoints.size() > 0) {
    // Restructure alpha and rhs b
    alpha.remove(deletedGridPoints);
    bNumP.remove(deletedGridPoints);
    bDenomP.remove(deletedGridPoints);
    bNumQ.remove(deletedGridPoints);
    bDenomQ.remove(deletedGridPoints);
  }
  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // Refinement, expand alpha and rhs b
  if (newNoPoints > oldNoPoints) {
    alpha.resizeZero(newNoPoints);
    bNumP.resizeZero(newNoPoints);
    bDenomP.resizeZero(newNoPoints);
    bNumQ.resizeZero(newNoPoints);
    bDenomQ.resizeZero(newNoPoints);
  }

  return true;
}

void ModelFittingDensityDifferenceEstimationCG::update(Dataset& newDatasetP, Dataset& newDatasetQ) {
  dataset = &newDatasetP;
  extraDataset = &newDatasetQ;
  update(newDatasetP.getData(), newDatasetQ.getData());
}

base::OperationMatrix* ModelFittingDensityDifferenceEstimationCG::computeRegularizationMatrix(
    base::Grid& grid) {
  base::OperationMatrix* C;
  auto& regularizationConfig = this->config->getRegularizationConfig();
  if (regularizationConfig.type_ == datadriven::RegularizationType::Identity) {
    C = op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.type_ == datadriven::RegularizationType::Laplace) {
    C = op_factory::createOperationLaplace(grid);
  } else {
    throw base::application_exception(
        "ModelFittingDensityDifferenceEstimationCG : unsupported regularization type");
  }
  return C;
}

void ModelFittingDensityDifferenceEstimationCG::update(DataMatrix& newDatasetP,
                                                       DataMatrix& newDatasetQ) {
  if (grid == nullptr) {
    // Initial fitting of datasets
    fit(newDatasetP, newDatasetQ);
  } else {
    // Compute regularization operator and retrieve the regularization strength
    auto& regularizationConfig = this->config->getRegularizationConfig();
    auto C = computeRegularizationMatrix(*grid);

    // Calculate the update for the rhs
    // Online procedure: beta is a forgetRate
    //    1 = forget all past batches
    //    0 = equal weighting
    DataVector rhsUpdate(grid->getSize());
    DataVector rhsUpdateExtra(grid->getSize());
    datadriven::DensityDifferenceSystemMatrix SMatrix(*grid, newDatasetP, newDatasetQ, C,
                                                      regularizationConfig.lambda_);
    SMatrix.computeUnweightedRhs(rhsUpdate, rhsUpdateExtra);
    double numInstancesP = static_cast<double>(newDatasetP.getNrows());
    double numInstancesQ = static_cast<double>(newDatasetQ.getNrows());
    double beta = this->config->getLearnerConfig().learningRate_;  // forgetRate
    // Update numerators
    bNumP.mult(1. - beta);
    bNumP.add(rhsUpdate);
    bNumQ.mult(1. - beta);
    bNumQ.add(rhsUpdateExtra);
    // Update denominators
    rhsUpdate.setAll(numInstancesP);
    bDenomP.mult(1. - beta);
    bDenomP.add(rhsUpdate);
    rhsUpdateExtra.setAll(numInstancesQ);
    bDenomQ.mult(1. - beta);
    bDenomQ.add(rhsUpdateExtra);
    // Compute the final rhs
    rhsUpdate = bNumP;
    rhsUpdate.componentwise_div(bDenomP);
    rhsUpdateExtra = bNumQ;
    rhsUpdateExtra.componentwise_div(bDenomQ);
    rhsUpdate.sub(rhsUpdateExtra);

    // Solve the system
    auto& solverConfig = this->config->getSolverRefineConfig();
    solver::ConjugateGradients cgSolver(solverConfig.maxIterations_, solverConfig.eps_);
    cgSolver.solve(SMatrix, alpha, rhsUpdate, true, solverConfig.verbose_, solverConfig.threshold_);
  }
}

bool ModelFittingDensityDifferenceEstimationCG::isRefinable() { return true; }

void ModelFittingDensityDifferenceEstimationCG::reset() {
  // Clear model
  grid.reset();
  refinementsPerformed = 0;
}

}  // namespace datadriven
}  // namespace sgpp
