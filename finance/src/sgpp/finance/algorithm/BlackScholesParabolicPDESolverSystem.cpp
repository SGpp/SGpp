// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/finance/operation/FinanceOpFactory.hpp>

#include <sgpp/globaldef.hpp>
#include <cmath>
#include <string>
#include <vector>

namespace SGPP {
namespace finance {

BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem(
    SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& mu,
    SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho, float_t r, float_t TimestepSize,
    std::string OperationMode, float_t dStrike, std::string option_type, bool bLogTransform,
    bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode, int numCoarsenPoints,
    float_t refineThreshold, std::string refineMode,
    SGPP::base::GridIndex::level_type refineMaxLevel) {
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;

  this->alpha_complete_old = new SGPP::base::DataVector(*this->alpha_complete);
  this->alpha_complete_tmp = new SGPP::base::DataVector(*this->alpha_complete);
  this->oldGridStorage = new SGPP::base::GridStorage(this->BoundGrid->getStorage());
  this->secondGridStorage = new SGPP::base::GridStorage(this->BoundGrid->getStorage());

  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->TimestepSize_old = TimestepSize;
  this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->r = r;
  this->mus = &mu;
  this->sigmas = &sigma;
  this->rhos = &rho;
  this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();

  // throw exception if grid dimensions not equal algorithmic dimensions
  if (this->BSalgoDims.size() > this->BoundGrid->getStorage().dim()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystemn : Number of "
        "algorithmic dimensions higher than the number of grid's dimensions.");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and
  // sigma)
  if (this->BoundGrid->getStorage().dim() != this->mus->getSize() ||
      this->BoundGrid->getStorage().dim() != this->sigmas->getSize()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Dimension of "
        "mu and sigma parameters don't match the grid's dimensions!");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (rho)
  if (this->BoundGrid->getStorage().dim() != this->rhos->getNrows() ||
      this->BoundGrid->getStorage().dim() != this->rhos->getNcols()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Row or col "
        "of rho parameter don't match the grid's dimensions!");
  }

  // test if all algorithmic dimensions are inside the grid's dimensions
  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    if (this->BSalgoDims[i] >= this->BoundGrid->getStorage().dim()) {
      throw SGPP::base::algorithm_exception(
          "BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Minimum "
          "one algorithmic dimension is not inside the grid's dimensions!");
    }
  }

  // test if there are float_t algorithmic dimensions
  std::vector<size_t> tempAlgoDims(this->BSalgoDims);

  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    size_t dimCount = 0;

    for (size_t j = 0; j < tempAlgoDims.size(); j++) {
      if (this->BSalgoDims[i] == tempAlgoDims[j]) {
        dimCount++;
      }
    }

    if (dimCount > 1) {
      throw SGPP::base::algorithm_exception(
          "BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : There is "
          "minimum one float_td algorithmic dimension!");
    }
  }

  // build the coefficient vectors for the operations
  this->gammaCoef = new SGPP::base::DataMatrix(this->BSalgoDims.size(), this->BSalgoDims.size());
  this->deltaCoef = new SGPP::base::DataVector(this->BSalgoDims.size());

  if (bLogTransform == false) {
    buildDeltaCoefficients();
    buildGammaCoefficients();

    // Create needed operations, on boundary grid
    this->OpDeltaBound =
        SGPP::op_factory::createOperationDelta(*this->BoundGrid, *this->deltaCoef).release();
    this->OpGammaBound =
        SGPP::op_factory::createOperationGamma(*this->BoundGrid, *this->gammaCoef).release();
  } else {
    // create needed operations that are different in case of a log-transformed Black-Scholoes
    // equation
    buildDeltaCoefficientsLogTransform();
    buildGammaCoefficientsLogTransform();

    // operations on boundary grid
    this->OpDeltaBound =
        SGPP::op_factory::createOperationDeltaLog(*this->BoundGrid, *this->deltaCoef).release();
    this->OpGammaBound =
        SGPP::op_factory::createOperationGammaLog(*this->BoundGrid, *this->gammaCoef).release();
  }

  this->OpLTwoBound = SGPP::op_factory::createOperationLTwoDotProduct(*this->BoundGrid).release();

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

  // init option type and strike
  this->dStrike = dStrike;
  this->option_type = option_type;

  // save coordinate transformations
  this->b_log_transform = bLogTransform;
}

BlackScholesParabolicPDESolverSystem::~BlackScholesParabolicPDESolverSystem() {
  delete this->OpDeltaBound;
  delete this->OpGammaBound;
  delete this->OpLTwoBound;
  delete this->gammaCoef;
  delete this->deltaCoef;
  delete this->BoundaryUpdate;

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  delete this->alpha_complete_old;
  delete this->alpha_complete_tmp;
  delete this->oldGridStorage;
  delete this->secondGridStorage;
}

void BlackScholesParabolicPDESolverSystem::applyLOperator(SGPP::base::DataVector& alpha,
                                                          SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the riskfree rate
  if (this->r != 0.0) {
    this->OpLTwoBound->mult(alpha, temp);
    result.axpy((-1.0) * this->r, temp);
  }

  // Apply the delta method
  this->OpDeltaBound->mult(alpha, temp);
  result.add(temp);

  // Apply the gamma method
  this->OpGammaBound->mult(alpha, temp);
  result.sub(temp);
}

void BlackScholesParabolicPDESolverSystem::applyMassMatrix(SGPP::base::DataVector& alpha,
                                                           SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoBound->mult(alpha, temp);

  result.add(temp);
}

void BlackScholesParabolicPDESolverSystem::finishTimestep() {
#ifndef NOBOUNDARYDISCOUNT

  // Adjust the boundaries with the riskfree rate
  if (this->r != 0.0) {
    if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas") {
      this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete,
                                             exp(((-1.0) * (this->r * this->TimestepSize))));
    }
  }

#endif

  // add number of Gridpoints
  this->numSumGridpointsInner += 0;
  this->numSumGridpointsComplete += this->BoundGrid->getSize();
}

void BlackScholesParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {
  if (this->useCoarsen == true && isLastTimestep == false) {
    ///////////////////////////////////////////////////
    // Start integrated refinement & coarsening
    ///////////////////////////////////////////////////

    size_t originalGridSize = this->BoundGrid->getStorage().size();

    // Coarsen the grid
    SGPP::base::GridGenerator& myGenerator = this->BoundGrid->getGenerator();

    // std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
    // std::cout << "Grid Size: " << originalGridSize << std::endl;

    if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numRefines = myGenerator.getNumberOfRefinablePoints();
      SGPP::base::SurplusRefinementFunctor myRefineFunc(
          *this->alpha_complete, numRefines, this->refineThreshold);

      if (this->refineMode == "maxLevel") {
        myGenerator.refineMaxLevel(myRefineFunc, this->refineMaxLevel);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage().size());
      }

      if (this->refineMode == "classic") {
        myGenerator.refine(myRefineFunc);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage().size());
      }
    }

    if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numCoarsen = myGenerator.getNumberOfRemovablePoints();
      SGPP::base::SurplusCoarseningFunctor myCoarsenFunctor(
          *this->alpha_complete, numCoarsen, this->coarsenThreshold);
      myGenerator.coarsenNFirstOnly(myCoarsenFunctor, *this->alpha_complete, originalGridSize);
    }

    ///////////////////////////////////////////////////
    // End integrated refinement & coarsening
    ///////////////////////////////////////////////////
  }
}

void BlackScholesParabolicPDESolverSystem::startTimestep() {
#ifndef NOBOUNDARYDISCOUNT

  // Adjust the boundaries with the riskfree rate
  if (this->r != 0.0) {
    if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
      this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete,
                                             exp(((-1.0) * (this->r * this->TimestepSize))));
    }
  }

#endif
}

void BlackScholesParabolicPDESolverSystem::buildGammaCoefficients() {
  size_t dim = this->BSalgoDims.size();

  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < dim; j++) {
      // handle diagonal
      if (i == j) {
        this->gammaCoef->set(i, j,
                             0.5 * ((this->sigmas->get(this->BSalgoDims[i]) *
                                     this->sigmas->get(this->BSalgoDims[j])) *
                                    this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
      } else {
        this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i]) *
                                     this->sigmas->get(this->BSalgoDims[j])) *
                                    this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
      }
    }
  }
}

void BlackScholesParabolicPDESolverSystem::buildDeltaCoefficients() {
  size_t dim = this->BSalgoDims.size();
  float_t covar_sum = 0.0;

  for (size_t i = 0; i < dim; i++) {
    covar_sum = 0.0;

    for (size_t j = 0; j < dim; j++) {
      // handle diagonal
      if (i == j) {
        covar_sum +=
            ((this->sigmas->get(this->BSalgoDims[i]) * this->sigmas->get(this->BSalgoDims[j])) *
             this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j]));
      } else {
        covar_sum +=
            (0.5 *
             ((this->sigmas->get(this->BSalgoDims[i]) * this->sigmas->get(this->BSalgoDims[j])) *
              this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
      }
    }

    this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i]) - covar_sum);
  }
}

void BlackScholesParabolicPDESolverSystem::buildGammaCoefficientsLogTransform() {
  size_t dim = this->BSalgoDims.size();

  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < dim; j++) {
      // handle diagonal
      if (i == j) {
        this->gammaCoef->set(i, j,
                             0.5 * ((this->sigmas->get(this->BSalgoDims[i]) *
                                     this->sigmas->get(this->BSalgoDims[j])) *
                                    this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
      } else {
        this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i]) *
                                     this->sigmas->get(this->BSalgoDims[j])) *
                                    this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
      }
    }
  }
}

void BlackScholesParabolicPDESolverSystem::buildDeltaCoefficientsLogTransform() {
  size_t dim = this->BSalgoDims.size();

  for (size_t i = 0; i < dim; i++) {
    this->deltaCoef->set(
        i, this->mus->get(this->BSalgoDims[i]) - (0.5 * (this->sigmas->get(this->BSalgoDims[i]) *
                                                         this->sigmas->get(this->BSalgoDims[i]))));
  }
}
}  // namespace finance
}  // namespace SGPP
