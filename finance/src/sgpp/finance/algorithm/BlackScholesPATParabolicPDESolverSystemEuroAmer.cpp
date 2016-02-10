// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesPATParabolicPDESolverSystemEuroAmer.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

BlackScholesPATParabolicPDESolverSystemEuroAmer::BlackScholesPATParabolicPDESolverSystemEuroAmer(
  SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha,
  SGPP::base::DataVector& lambda,
  SGPP::base::DataMatrix& eigenvecs, SGPP::base::DataVector& mu_hat,
  float_t TimestepSize, std::string OperationMode,
  float_t dStrike, std::string option_type, float_t r,
  bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode,
  int numCoarsenPoints, float_t refineThreshold, std::string refineMode,
  SGPP::base::GridIndex::level_type refineMaxLevel) {
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;

  this->alpha_complete_old = new SGPP::base::DataVector(*this->alpha_complete);
  this->alpha_complete_tmp = new SGPP::base::DataVector(*this->alpha_complete);
  this->oldGridStorage = new SGPP::base::GridStorage(*
      (this->BoundGrid)->getStorage());
  this->secondGridStorage = new SGPP::base::GridStorage(*
      (this->BoundGrid)->getStorage());

  this->InnerGrid = NULL;
  this->alpha_inner = NULL;
  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->TimestepSize_old = TimestepSize;
  this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(
    SparseGrid.getStorage());
  this->GridConverter = new SGPP::base::DirichletGridConverter();

  // set Eigenvalues, Eigenvector of covariance matrix and mu_hat
  this->lambda = new SGPP::base::DataVector(lambda);
  this->eigenvecs = new SGPP::base::DataMatrix(eigenvecs);
  this->mu_hat = new SGPP::base::DataVector(mu_hat);

  this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();
  this->nExecTimesteps = 0;

  // throw exception if grid dimensions not equal algorithmic dimensions
  if (this->BSalgoDims.size() > this->BoundGrid->getStorage()->dim()) {
    throw SGPP::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystemn : Number of algorithmic dimensions higher than the number of grid's dimensions.");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and sigma)
  if (this->BoundGrid->getStorage()->dim() != this->lambda->getSize()) {
    throw SGPP::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : Dimension of mu and sigma parameters don't match the grid's dimensions!");
  }

  // test if all algorithmic dimensions are inside the grid's dimensions
  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    if (this->BSalgoDims[i] >= this->BoundGrid->getStorage()->dim()) {
      throw SGPP::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : Minimum one algorithmic dimension is not inside the grid's dimensions!");
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
      throw SGPP::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : There is minimum one float_td algorithmic dimension!");
    }
  }

  // create the inner grid
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid,
      *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

  // Create operations
  this->OpLaplaceInner = SGPP::op_factory::createOperationLaplace(
                           *this->InnerGrid, *this->lambda);
  this->OpLaplaceBound = SGPP::op_factory::createOperationLaplace(
                           *this->BoundGrid, *this->lambda);

  this->OpLTwoInner = SGPP::op_factory::createOperationLTwoDotProduct(
                        *this->InnerGrid);
  this->OpLTwoBound = SGPP::op_factory::createOperationLTwoDotProduct(
                        *this->BoundGrid);

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
  this->r = r;
}

BlackScholesPATParabolicPDESolverSystemEuroAmer::~BlackScholesPATParabolicPDESolverSystemEuroAmer() {
  delete this->OpLaplaceBound;
  delete this->OpLTwoBound;
  delete this->OpLaplaceInner;
  delete this->OpLTwoInner;
  delete this->BoundaryUpdate;
  delete this->GridConverter;

  if (this->InnerGrid != NULL) {
    delete this->InnerGrid;
  }

  if (this->alpha_inner != NULL) {
    delete this->alpha_inner;
  }

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  delete this->alpha_complete_old;
  delete this->alpha_complete_tmp;
  delete this->lambda;
  delete this->eigenvecs;
  delete this->mu_hat;
}

void BlackScholesPATParabolicPDESolverSystemEuroAmer::applyLOperatorComplete(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the Laplace operator
  this->OpLaplaceBound->mult(alpha, temp);
  result.axpy(-0.5, temp);
}

void BlackScholesPATParabolicPDESolverSystemEuroAmer::applyLOperatorInner(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the Laplace operator
  this->OpLaplaceInner->mult(alpha, temp);
  result.axpy(-0.5, temp);
}

void BlackScholesPATParabolicPDESolverSystemEuroAmer::applyMassMatrixComplete(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoBound->mult(alpha, temp);

  result.add(temp);
}

void BlackScholesPATParabolicPDESolverSystemEuroAmer::applyMassMatrixInner(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoInner->mult(alpha, temp);

  result.add(temp);
}


void BlackScholesPATParabolicPDESolverSystemEuroAmer::finishTimestep() {
  // Replace the inner coefficients on the boundary grid
  this->GridConverter->updateBoundaryCoefs(*this->alpha_complete,
      *this->alpha_inner);

  // check if we are doing an American put -> handle early exercise
  if (this->option_type == "std_amer_put") {
    float_t current_time = static_cast<float_t>(this->nExecTimesteps) *
                           this->TimestepSize;

    SGPP::base::OperationHierarchisation* myHierarchisation =
      SGPP::op_factory::createOperationHierarchisation(*this->BoundGrid);
    myHierarchisation->doDehierarchisation(*this->alpha_complete);
    size_t dim = this->BoundGrid->getStorage()->dim();
    SGPP::base::BoundingBox* myBB = new SGPP::base::BoundingBox(*
        (this->BoundGrid->getBoundingBox()));

    float_t* coords_val = new float_t[dim];

    for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++) {
      std::vector<float_t> eval_point_coord;
      std::string coords = this->BoundGrid->getStorage()->get(i)->getCoordsStringBB(
                             *myBB);
      std::stringstream coordsStream(coords);

      float_t tmp;

      // read coordinates
      for (size_t j = 0; j < dim; j++) {
        coordsStream >> tmp;

        coords_val[j] = tmp;
      }

      tmp = 0.0;

      for (size_t j = 0; j < dim; j++) {
        float_t inner_tmp = 0.0;

        for (size_t l = 0; l < dim; l++) {
          inner_tmp += this->eigenvecs->get(j,
                                            l) * (coords_val[l] - (current_time * this->mu_hat->get(l)));
        }

        tmp += exp(inner_tmp);
      }

      float_t payoff = std::max<float_t>(this->dStrike - (tmp / static_cast<float_t>
                                         (dim)), 0.0);
      float_t discounted_value = ((*this->alpha_complete)[i]) * exp(((-1.0) *
                                 (this->r * this->TimestepSize)));

      (*this->alpha_complete)[i] = std::max<float_t>(payoff, discounted_value);
    }

    delete[] coords_val;

    myHierarchisation->doHierarchisation(*this->alpha_complete);
    delete myHierarchisation;
    delete myBB;
  }

}
void BlackScholesPATParabolicPDESolverSystemEuroAmer::coarsenAndRefine(
  bool isLastTimestep) {
  // add number of Gridpoints
  this->numSumGridpointsInner += this->InnerGrid->getSize();
  this->numSumGridpointsComplete += this->BoundGrid->getSize();

  if (this->useCoarsen == true && isLastTimestep == false) {
    ///////////////////////////////////////////////////
    // Start integrated refinement & coarsening
    ///////////////////////////////////////////////////

    size_t originalGridSize = this->BoundGrid->getStorage()->size();

    // Coarsen the grid
    SGPP::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

    //std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
    //std::cout << "Grid Size: " << originalGridSize << std::endl;

    if (this->adaptSolveMode == "refine"
        || this->adaptSolveMode == "coarsenNrefine") {
      size_t numRefines = myGenerator->getNumberOfRefinablePoints();
      SGPP::base::SurplusRefinementFunctor* myRefineFunc = new
      SGPP::base::SurplusRefinementFunctor(this->alpha_complete, numRefines,
                                           this->refineThreshold);

      if (this->refineMode == "maxLevel") {
        myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
      }

      if (this->refineMode == "classic") {
        myGenerator->refine(myRefineFunc);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
      }

      delete myRefineFunc;
    }

    if (this->adaptSolveMode == "coarsen"
        || this->adaptSolveMode == "coarsenNrefine") {
      size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
      SGPP::base::SurplusCoarseningFunctor* myCoarsenFunctor = new
      SGPP::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen,
                                           this->coarsenThreshold);
      myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete,
                                     originalGridSize);
      delete myCoarsenFunctor;
    }

    delete myGenerator;

    ///////////////////////////////////////////////////
    // End integrated refinement & coarsening
    ///////////////////////////////////////////////////

    // rebuild the inner grid + coefficients
    this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid,
        *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
  }

}

void BlackScholesPATParabolicPDESolverSystemEuroAmer::startTimestep() {
}

}

}