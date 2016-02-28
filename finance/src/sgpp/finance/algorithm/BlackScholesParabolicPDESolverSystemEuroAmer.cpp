// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/finance/operation/FinanceOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

namespace SGPP {
namespace finance {

BlackScholesParabolicPDESolverSystemEuroAmer::BlackScholesParabolicPDESolverSystemEuroAmer(
    SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& mu,
    SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho, float_t r, float_t TimestepSize,
    std::string OperationMode, float_t dStrike, std::string option_type, bool bLogTransform,
    bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode, int numCoarsenPoints,
    float_t refineThreshold, std::string refineMode, SGPP::base::GridIndex::level_type xLevel) {
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;

  this->alpha_complete_old = new SGPP::base::DataVector(*this->alpha_complete);
  this->alpha_complete_tmp = new SGPP::base::DataVector(*this->alpha_complete);
  this->oldGridStorage = new SGPP::base::GridStorage(this->BoundGrid->getStorage());
  this->secondGridStorage = new SGPP::base::GridStorage(this->BoundGrid->getStorage());

  this->InnerGrid = NULL;
  this->alpha_inner = NULL;
  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->TimestepSize_old = TimestepSize;
  this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->GridConverter = new SGPP::base::DirichletGridConverter();
  this->r = r;
  this->mus = &mu;
  this->sigmas = &sigma;
  this->rhos = &rho;
  this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();
  this->nExecTimesteps = 0;

  // throw exception if grid dimensions not equal algorithmic dimensions
  if (this->BSalgoDims.size() != this->BoundGrid->getDimension()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystemEuropean::"
        "BlackScholesParabolicPDESolverSystemEuropean : Number of algorithmic dimensions is not "
        "equal to the number of grid's dimensions.");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and
  // sigma)
  if (this->BoundGrid->getDimension() != this->mus->getSize() ||
      this->BoundGrid->getDimension() != this->sigmas->getSize()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystemEuropean::"
        "BlackScholesParabolicPDESolverSystemEuropean : Dimension of mu and sigma parameters don't "
        "match the grid's dimensions!");
  }

  // test if number of dimensions in the coefficients match the numbers of grid dimensions (rho)
  if (this->BoundGrid->getDimension() != this->rhos->getNrows() ||
      this->BoundGrid->getDimension() != this->rhos->getNcols()) {
    throw SGPP::base::algorithm_exception(
        "BlackScholesParabolicPDESolverSystemEuropean::"
        "BlackScholesParabolicPDESolverSystemEuropean : Row or col of rho parameter don't match "
        "the grid's dimensions!");
  }

  // test if all algorithmic dimensions are inside the grid's dimensions
  for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
    if (this->BSalgoDims[i] >= this->BoundGrid->getDimension()) {
      throw SGPP::base::algorithm_exception(
          "BlackScholesParabolicPDESolverSystemEuropean::"
          "BlackScholesParabolicPDESolverSystemEuropean : Minimum one algorithmic dimension is not "
          "inside the grid's dimensions!");
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
          "BlackScholesParabolicPDESolverSystemEuropean::"
          "BlackScholesParabolicPDESolverSystemEuropean : There is minimum one float_td "
          "algorithmic dimension!");
    }
  }

  // build the coefficient vectors for the operations
  this->gammaCoef = new SGPP::base::DataMatrix(this->BSalgoDims.size(), this->BSalgoDims.size());
  this->deltaCoef = new SGPP::base::DataVector(this->BSalgoDims.size());

  // create the inner grid
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                               &this->InnerGrid, &this->alpha_inner);

  if (bLogTransform == false) {
    buildDeltaCoefficients();
    buildGammaCoefficients();

    // Create needed operations, on inner grid
    this->OpDeltaInner =
        SGPP::op_factory::createOperationDelta(*this->InnerGrid, *this->deltaCoef).release();
    this->OpGammaInner =
        SGPP::op_factory::createOperationGamma(*this->InnerGrid, *this->gammaCoef).release();
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
    // operations on inner grid
    this->OpDeltaInner =
        SGPP::op_factory::createOperationDeltaLog(*this->InnerGrid, *this->deltaCoef).release();
    this->OpGammaInner =
        SGPP::op_factory::createOperationGammaLog(*this->InnerGrid, *this->gammaCoef).release();
  }

  // Create operations, independent bLogTransform
  this->OpLTwoInner = SGPP::op_factory::createOperationLTwoDotProduct(*this->InnerGrid).release();
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

#ifdef HEDGE
  SGPP::base::BoundingBox* grid_bb = this->BoundGrid->getBoundingBox();
  SGPP::base::DimensionBoundary* myBoundaries =
      new SGPP::base::DimensionBoundary[grid_bb->getDimensions()];

  for (size_t d = 0; d < grid_bb->getDimensions(); d++) {
    if (bLogTransform == true) {
      float_t interval_width =
          exp(grid_bb->getBoundary(d).rightBoundary) - exp(grid_bb->getBoundary(d).leftBoundary);
      float_t hedge_offset = (interval_width - (interval_width * HEDGE_WIDTH_PERCENT)) / 2.0;
      myBoundaries[d].leftBoundary = exp(grid_bb->getBoundary(d).leftBoundary) + hedge_offset;
      myBoundaries[d].rightBoundary = exp(grid_bb->getBoundary(d).rightBoundary) - hedge_offset;
    } else {
      float_t hedge_offset =
          (grid_bb->getIntervalWidth(d) - (grid_bb->getIntervalWidth(d) * HEDGE_WIDTH_PERCENT)) /
          2.0;
      myBoundaries[d].leftBoundary = grid_bb->getBoundary(d).leftBoundary + hedge_offset;
      myBoundaries[d].rightBoundary = grid_bb->getBoundary(d).rightBoundary - hedge_offset;
    }

    myBoundaries[d].bDirichletLeft = true;
    myBoundaries[d].bDirichletRight = true;
  }

  SGPP::base::BoundingBox* myHedgeBB =
      new SGPP::base::BoundingBox(grid_bb->getDimensions(), myBoundaries);
  // hedging
  myHedge = new SGPP::finance::Hedging(*myHedgeBB, HEDGE_POINTS_PER_DIM, HEDGE_EPS, bLogTransform);

  delete myHedgeBB;
  delete[] myBoundaries;
#endif
}

BlackScholesParabolicPDESolverSystemEuroAmer::~BlackScholesParabolicPDESolverSystemEuroAmer() {
  delete this->OpDeltaBound;
  delete this->OpGammaBound;
  delete this->OpLTwoBound;
  delete this->OpDeltaInner;
  delete this->OpGammaInner;
  delete this->OpLTwoInner;
  delete this->gammaCoef;
  delete this->deltaCoef;
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

#ifdef HEDGE
  delete myHedge;
#endif
}

void BlackScholesParabolicPDESolverSystemEuroAmer::applyLOperatorComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
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

void BlackScholesParabolicPDESolverSystemEuroAmer::applyLOperatorInner(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the riskfree rate
  if (this->r != 0.0) {
    this->OpLTwoInner->mult(alpha, temp);
    result.axpy((-1.0) * this->r, temp);
  }

  // Apply the delta method
  this->OpDeltaInner->mult(alpha, temp);
  result.add(temp);

  // Apply the gamma method
  this->OpGammaInner->mult(alpha, temp);
  result.sub(temp);
}

void BlackScholesParabolicPDESolverSystemEuroAmer::applyMassMatrixComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoBound->mult(alpha, temp);

  result.add(temp);
}

void BlackScholesParabolicPDESolverSystemEuroAmer::applyMassMatrixInner(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoInner->mult(alpha, temp);

  result.add(temp);
}

void BlackScholesParabolicPDESolverSystemEuroAmer::finishTimestep() {
  this->nExecTimesteps++;

  // Replace the inner coefficients on the boundary grid
  this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

#ifndef NOBOUNDARYDISCOUNT

  // Adjust the boundaries with the riskfree rate
  if (this->r != 0.0) {
    if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas") {
      this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete,
                                             exp(((-1.0) * (this->r * this->TimestepSize))));
    }
  }

#endif

  // check if we are doing an American put -> handle early exercise
  if (this->option_type == "std_amer_put") {
    std::unique_ptr<SGPP::base::OperationHierarchisation> myHierarchisation(
        SGPP::op_factory::createOperationHierarchisation(*this->BoundGrid));
    myHierarchisation->doDehierarchisation(*this->alpha_complete);
    size_t dim = this->BoundGrid->getDimension();
    SGPP::base::BoundingBox* myBB =
        new SGPP::base::BoundingBox(this->BoundGrid->getBoundingBox());

    float_t* dblFuncValues = new float_t[dim];

    for (size_t i = 0; i < this->BoundGrid->getSize(); i++) {
      std::string coords = this->BoundGrid->getStorage().get(i)->getCoordsStringBB(*myBB);
      std::stringstream coordsStream(coords);

      float_t tmp;

      // read coordinates
      for (size_t j = 0; j < dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      tmp = 0.0;

      if (this->b_log_transform == true) {
        for (size_t j = 0; j < dim; j++) {
          tmp += exp(dblFuncValues[j]);
        }
      } else {
        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }
      }

      //      (*this->alpha_complete)[i] = std::max<float_t>((*this->alpha_complete)[i],
      //      (std::max<float_t>(this->dStrike-((tmp/static_cast<float_t>(dim))),
      //      0.0))*exp(((-1.0)*(this->r*
      //      static_cast<float_t>(this->nExecTimesteps)*this->TimestepSize))));
      (*this->alpha_complete)[i] = std::max<float_t>(
          (*this->alpha_complete)[i],
          (std::max<float_t>(this->dStrike - ((tmp / static_cast<float_t>(dim))), 0.0)));
    }

    delete[] dblFuncValues;

    myHierarchisation->doHierarchisation(*this->alpha_complete);
    delete myBB;
  }
}

void BlackScholesParabolicPDESolverSystemEuroAmer::coarsenAndRefine(bool isLastTimestep) {
  // add number of Gridpoints
  this->numSumGridpointsInner += this->InnerGrid->getSize();
  this->numSumGridpointsComplete += this->BoundGrid->getSize();

  if (this->useCoarsen == true && isLastTimestep == false) {
    ///////////////////////////////////////////////////
    // Start integrated refinement & coarsening
    ///////////////////////////////////////////////////

    size_t originalGridSize = this->BoundGrid->getSize();

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
        this->alpha_complete->resizeZero(this->BoundGrid->getSize());
      }

      if (this->refineMode == "classic") {
        myGenerator.refine(myRefineFunc);
        this->alpha_complete->resizeZero(this->BoundGrid->getSize());
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

    // rebuild the inner grid + coefficients
    this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                                   &this->InnerGrid, &this->alpha_inner);
  }

#ifdef HEDGE
  std::stringstream filename_ext;
  filename_ext << ((this->nExecTimesteps) * this->TimestepSize);

  myHedge->calc_hedging(*this->BoundGrid, *this->alpha_complete, filename_ext.str());
#endif
}

void BlackScholesParabolicPDESolverSystemEuroAmer::startTimestep() {
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

void BlackScholesParabolicPDESolverSystemEuroAmer::buildGammaCoefficients() {
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

void BlackScholesParabolicPDESolverSystemEuroAmer::buildDeltaCoefficients() {
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

void BlackScholesParabolicPDESolverSystemEuroAmer::buildGammaCoefficientsLogTransform() {
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

void BlackScholesParabolicPDESolverSystemEuroAmer::buildDeltaCoefficientsLogTransform() {
  size_t dim = this->BSalgoDims.size();

  for (size_t i = 0; i < dim; i++) {
    this->deltaCoef->set(
        i, this->mus->get(this->BSalgoDims[i]) - (0.5 * (this->sigmas->get(this->BSalgoDims[i]) *
                                                         this->sigmas->get(this->BSalgoDims[i]))));
  }
}
}  // namespace finance
}  // namespace SGPP
