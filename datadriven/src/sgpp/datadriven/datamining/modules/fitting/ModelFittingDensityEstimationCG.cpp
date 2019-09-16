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
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <string>
#include <vector>
#include <list>


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

ModelFittingDensityEstimationCG::ModelFittingDensityEstimationCG(
    const FitterConfigurationDensityEstimation& config) : ModelFittingDensityEstimation{},
        bNum{0}, bDenom{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityEstimationCG::evaluate(const DataVector& sample) {
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
  return opEval->eval(alpha, sample);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityEstimationCG::evaluate(DataMatrix& samples, DataVector& results) {
  sgpp::op_factory::createOperationMultipleEval(*grid, samples)->eval(alpha, results);
}

void ModelFittingDensityEstimationCG::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationCG::fit(DataMatrix& newDataset) {
  // Reset grid
  reset();

  // Setup new grid
  auto& gridConfig = this->config->getGridConfig();
  auto& geometryConfig = this->config->getGeometryConfig();
  gridConfig.dim_ = newDataset.getNcols();
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig, geometryConfig)};
  // build surplus vector
  alpha = DataVector{grid->getSize()};

  // Initialize the right hand side (numerator and denominator)
  bNum = DataVector{grid->getSize()};
  bDenom = DataVector{grid->getSize()};
  bNum.setAll(0.0);
  bDenom.setAll(0.0);

  // Now that everythin is setup with zero data, simply call the update method
  update(newDataset);
}


bool ModelFittingDensityEstimationCG::refine(size_t newNoPoints,
    std::list<size_t> *deletedGridPoints) {
  // Coarsening, remove idx from alpha
  if (deletedGridPoints != nullptr && deletedGridPoints->size() > 0) {
    // Restructure alpha and rhs b
    std::vector<size_t> idxToDelete{std::begin(*deletedGridPoints), std::end(*deletedGridPoints)};
    alpha.remove(idxToDelete);
    bNum.remove(idxToDelete);
    bDenom.remove(idxToDelete);
  }
  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // Refinement, expand alpha and rhs b
  if (newNoPoints > oldNoPoints) {
    alpha.resizeZero(newNoPoints);
    bNum.resizeZero(newNoPoints);
    bDenom.resizeZero(newNoPoints);
  }

  return true;
}

void ModelFittingDensityEstimationCG::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

base::OperationMatrix* ModelFittingDensityEstimationCG::computeRegularizationMatrix(
    base::Grid& grid) {
  base::OperationMatrix* C;
  auto& regularizationConfig = this->config->getRegularizationConfig();
  if (regularizationConfig.type_ == datadriven::RegularizationType::Identity) {
    C = op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.type_ == datadriven::RegularizationType::Laplace) {
    C = op_factory::createOperationLaplace(grid);
  } else {
    throw base::application_exception(
        "ModelFittingDensityEstimationCG : unsupported regularization type");
  }
  return C;
}

void ModelFittingDensityEstimationCG::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Compute regularization operator and retrieve the regularization strength
    auto& regularizationConfig = this->config->getRegularizationConfig();
    auto C = computeRegularizationMatrix(*grid);

    // Calculate the update for the rhs
    DataVector rhsUpdate(grid->getSize());
    datadriven::DensitySystemMatrix SMatrix(*grid, newDataset, C, regularizationConfig.lambda_);
    SMatrix.generateb(rhsUpdate);
    double numInstances = static_cast<double>(newDataset.getNrows());
    // Rescale the rhs such that it is not normalized by the number of instances
    rhsUpdate.mult(static_cast<double>(numInstances));
    // Weigh the current right hand side with beta (decay)
    bNum.mult(this->config->getLearnerConfig().beta);

    bNum.add(rhsUpdate);
    // Update the denominator (dataset size) as well
    rhsUpdate.setAll(numInstances);
    bDenom.add(rhsUpdate);
    // Compute the current rhs
    rhsUpdate.setAll(0.0);
    rhsUpdate.add(bNum);
    rhsUpdate.componentwise_div(bDenom);

    // Solve the system
    auto& solverConfig = this->config->getSolverRefineConfig();
    solver::ConjugateGradients cgSolver(solverConfig.maxIterations_, solverConfig.eps_);
    cgSolver.solve(SMatrix, alpha, rhsUpdate, true, solverConfig.verbose_, solverConfig.threshold_);
  }
}

bool ModelFittingDensityEstimationCG::isRefinable() { return true; }

void ModelFittingDensityEstimationCG::reset() {
  // Clear model
  grid.reset();
  refinementsPerformed = 0;
}

}  // namespace datadriven
}  // namespace sgpp
