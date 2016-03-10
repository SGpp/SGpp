// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesPATParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace sgpp {

namespace finance {

BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, sgpp::base::DataVector& lambda,
    sgpp::base::DataMatrix& eigenvecs, sgpp::base::DataVector& mu_hat, double TimestepSize,
    std::string OperationMode, double dStrike, std::string option_type, bool useCoarsen,
    double coarsenThreshold, std::string adaptSolveMode, int numCoarsenPoints,
    double refineThreshold, std::string refineMode,
    sgpp::base::GridIndex::level_type refineMaxLevel) {
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;

  this->alpha_complete_old = new sgpp::base::DataVector(*this->alpha_complete);
  this->alpha_complete_tmp = new sgpp::base::DataVector(*this->alpha_complete);
  this->oldGridStorage = new sgpp::base::GridStorage(this->BoundGrid->getStorage());
  this->secondGridStorage = new sgpp::base::GridStorage(this->BoundGrid->getStorage());

  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->TimestepSize_old = TimestepSize;
  this->BoundaryUpdate = new sgpp::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();

  // set Eigenvalues, Eigenvector of covariance matrix and mu_hat
  this->lambda = new sgpp::base::DataVector(lambda);
  this->eigenvecs = new sgpp::base::DataMatrix(eigenvecs);
  this->mu_hat = new sgpp::base::DataVector(mu_hat);

  // throw exception if grid dimensions not equal algorithmic dimensions
  if (this->BSalgoDims.size() > this->BoundGrid->getDimension()) {
    throw sgpp::base::algorithm_exception(
        "BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystemn : "
        "Number of algorithmic dimensions higher than the number of grid's dimensions.");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and
  // sigma)
  if (this->BoundGrid->getDimension() != this->lambda->getSize()) {
    throw sgpp::base::algorithm_exception(
        "BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : "
        "Dimension of mu and sigma parameters don't match the grid's dimensions!");
  }

  // test if all algorithmic dimensions are inside the grid's dimensions
  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    if (this->BSalgoDims[i] >= this->BoundGrid->getDimension()) {
      throw sgpp::base::algorithm_exception(
          "BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : "
          "Minimum one algorithmic dimension is not inside the grid's dimensions!");
    }
  }

  // test if there are double algorithmic dimensions
  std::vector<size_t> tempAlgoDims(this->BSalgoDims);

  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    size_t dimCount = 0;

    for (size_t j = 0; j < tempAlgoDims.size(); j++) {
      if (this->BSalgoDims[i] == tempAlgoDims[j]) {
        dimCount++;
      }
    }

    if (dimCount > 1) {
      throw sgpp::base::algorithm_exception(
          "BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : "
          "There is minimum one doubled algorithmic dimension!");
    }
  }

  // operations on boundary grid
  this->OpLaplaceBound =
      sgpp::op_factory::createOperationLaplace(*this->BoundGrid, *this->lambda).release();
  this->OpLTwoBound =
      sgpp::op_factory::createOperationLTwoDotProduct(*this->BoundGrid).release();

  // right hand side if System
  this->rhs = NULL;

  // set coarsen settings
  this->useCoarsen = useCoarsen;
  this->coarsenThreshold = coarsenThreshold;
  this->refineThreshold = refineThreshold;
  this->adaptSolveMode = adaptSolveMode;
  this->numCoarsenPoints = numCoarsenPoints;
  this->refineMode = refineMode;
  this->refineMaxLevel = refineMaxLevel;

  // init Number of AverageGridPoins
  this->numSumGridpointsInner = 0;
  this->numSumGridpointsComplete = 0;
}

BlackScholesPATParabolicPDESolverSystem::~BlackScholesPATParabolicPDESolverSystem() {
  delete this->OpLaplaceBound;
  delete this->OpLTwoBound;
  delete this->BoundaryUpdate;

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  delete this->alpha_complete_old;
  delete this->alpha_complete_tmp;
  delete this->lambda;
  delete this->eigenvecs;
  delete this->mu_hat;
}

void BlackScholesPATParabolicPDESolverSystem::applyLOperator(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());
  result.setAll(0.0);

  // Apply the Laplace operator
  this->OpLaplaceBound->mult(alpha, temp);
  result.axpy(-0.5, temp);
}

void BlackScholesPATParabolicPDESolverSystem::applyMassMatrix(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());
  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoBound->mult(alpha, temp);

  result.add(temp);
}

void BlackScholesPATParabolicPDESolverSystem::finishTimestep() {}

void BlackScholesPATParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {
  // add number of Gridpoints
  this->numSumGridpointsInner += 0;
  this->numSumGridpointsComplete += this->BoundGrid->getSize();

  if (this->useCoarsen == true && isLastTimestep == false) {
    ///////////////////////////////////////////////////
    // Start integrated refinement & coarsening
    ///////////////////////////////////////////////////

    size_t originalGridSize = this->BoundGrid->getSize();

    // Coarsen the grid
    sgpp::base::GridGenerator& myGenerator = this->BoundGrid->getGenerator();

    // std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
    // std::cout << "Grid Size: " << originalGridSize << std::endl;

    if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numRefines = myGenerator.getNumberOfRefinablePoints();
      sgpp::base::SurplusRefinementFunctor myRefineFunc(
          *this->alpha_complete, numRefines, this->refineThreshold);

      if (this->refineMode == "maxLevel") {
        myGenerator.refineMaxLevel(myRefineFunc, this->refineMaxLevel);
        this->alpha_complete->resizeZero(this->BoundGrid->getSize());
      }

      if (this->refineMode == "classic") {
        myGenerator.refine(myRefineFunc);
        this->alpha_complete->resizeZero(this->BoundGrid->getSize());
      }
    }

    if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numCoarsen = myGenerator.getNumberOfRemovablePoints();
      sgpp::base::SurplusCoarseningFunctor myCoarsenFunctor(
          *this->alpha_complete, numCoarsen, this->coarsenThreshold);
      myGenerator.coarsenNFirstOnly(myCoarsenFunctor, *this->alpha_complete, originalGridSize);
    }

    ///////////////////////////////////////////////////
    // End integrated refinement & coarsening
    ///////////////////////////////////////////////////
  }
}

void BlackScholesPATParabolicPDESolverSystem::startTimestep() {}
}  // namespace finance
}  // namespace sgpp
