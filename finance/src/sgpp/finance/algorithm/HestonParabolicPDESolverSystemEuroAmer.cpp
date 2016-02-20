// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/HestonParabolicPDESolverSystemEuroAmer.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/finance/operation/FinanceOpFactory.hpp>

// #define NOBOUNDARYDISCOUNT

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace SGPP {
namespace finance {

HestonParabolicPDESolverSystemEuroAmer::HestonParabolicPDESolverSystemEuroAmer(
    SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& thetas,
    SGPP::base::DataVector& volvols, SGPP::base::DataVector& kappas, SGPP::base::DataMatrix& rho,
    float_t r, float_t TimestepSize, std::string OperationMode, float_t dStrike,
    std::string option_type, bool bLogTransform, bool useCoarsen, float_t coarsenThreshold,
    std::string adaptSolveMode, int numCoarsenPoints, float_t refineThreshold,
    std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel) {
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
  this->thetas = &thetas;
  this->volvols = &volvols;
  this->kappas = &kappas;
  this->hMatrix = &rho;
  this->HestonAlgoDims = this->BoundGrid->getAlgorithmicDimensions();
  this->nAssets = this->HestonAlgoDims.size() / 2;
  this->nExecTimesteps = 0;

  // throw an exception if we're not using Log coordinates and there is more than one asset
  if (this->nAssets > 1 && !bLogTransform) {
    throw SGPP::base::application_exception(
        "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : "
        "Cartesian coordinates are not supported for the basket case.");
  }

  // throw exception if algorithmic dimensions isn't even
  if ((this->HestonAlgoDims.size() % 2) != 0) {
    throw SGPP::base::algorithm_exception(
        "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Number "
        "of algorithmic dimensions is not even.");
  }

  // throw exception if grid dimensions not equal algorithmic dimensions
  if (this->HestonAlgoDims.size() != this->BoundGrid->getStorage().dim()) {
    throw SGPP::base::algorithm_exception(
        "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Number "
        "of algorithmic dimensions is not equal to the number of grid's dimensions.");
  }

  // Test if 2*dimC = dimG, where dimC is the number of dimensions in the coefficient vectors and
  // dimG is the number of grid dimensions
  if (this->BoundGrid->getStorage().dim() != (2 * this->thetas->getSize()) ||
      this->BoundGrid->getStorage().dim() != (2 * this->kappas->getSize()) ||
      this->BoundGrid->getStorage().dim() != (2 * this->volvols->getSize())) {
    throw SGPP::base::algorithm_exception(
        "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : "
        "Dimension of theta/volvol/kappa parameters != half of grid's dimensions!");
  }

  // Test if number of dimensions in the coefficients match the numbers of grid dimensions (hmatrix)
  if (this->BoundGrid->getStorage().dim() != this->hMatrix->getNrows() ||
      this->BoundGrid->getStorage().dim() != this->hMatrix->getNcols()) {
    throw SGPP::base::algorithm_exception(
        "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Row or "
        "col of hmatrix parameter don't match the grid's dimensions!");
  }

  // Test if all algorithmic dimensions are inside the grid's dimensions
  for (size_t i = 0; i < this->HestonAlgoDims.size(); i++) {
    if (this->HestonAlgoDims[i] >= this->BoundGrid->getStorage().dim()) {
      throw SGPP::base::algorithm_exception(
          "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : "
          "Minimum one algorithmic dimension is not inside the grid's dimensions!");
    }
  }

  // Test if there are float_t algorithmic dimensions
  std::vector<size_t> tempAlgoDims(this->HestonAlgoDims);

  for (size_t i = 0; i < this->HestonAlgoDims.size(); i++) {
    size_t dimCount = 0;

    for (size_t j = 0; j < tempAlgoDims.size(); j++) {
      if (this->HestonAlgoDims[i] == tempAlgoDims[j]) {
        dimCount++;
      }
    }

    if (dimCount > 1) {
      throw SGPP::base::algorithm_exception(
          "HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : There "
          "is minimum one float_td algorithmic dimension!");
    }
  }

  // Build the coefficient collections for the operations
  // Some of the coefficient collections are vectors, some are matrices and one is a tensor (for the
  // Up/Down Four Op Dim case)
  size_t coefficientVectorSize = this->HestonAlgoDims.size();
  this->dCoeff = new SGPP::base::DataVector(coefficientVectorSize);
  this->eCoeff = new SGPP::base::DataVector(coefficientVectorSize);
  this->fCoeff = new SGPP::base::DataVector(coefficientVectorSize);
  this->gCoeff = new SGPP::base::DataVector(coefficientVectorSize);
  this->zCoeff = new SGPP::base::DataVector(coefficientVectorSize);
  this->bCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  this->cCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  this->hCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  this->xCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  this->yCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  this->wCoeff = new SGPP::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
  create4dEqualDimSizeArray(coefficientVectorSize, &(this->kCoeff));

  // create the inner grid
  // Here we have to make sure that the alpha_inner doesn't come out with zeros at the other end.
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                               &this->InnerGrid, &this->alpha_inner);

  if (bLogTransform == false) {
    // Build coefficients
    buildXCoefficients();
    buildYCoefficients();
    buildWCoefficients();
    buildZCoefficients();
    buildGCoefficients();
    buildFCoefficients();
    buildDCoefficients();

    // Create operators
    this->OpXBound = SGPP::op_factory::createOperationHestonX(*this->BoundGrid, *this->xCoeff);
    this->OpXInner = SGPP::op_factory::createOperationHestonX(*this->InnerGrid, *this->xCoeff);
    this->OpYBound = SGPP::op_factory::createOperationHestonY(*this->BoundGrid, *this->yCoeff);
    this->OpYInner = SGPP::op_factory::createOperationHestonY(*this->InnerGrid, *this->yCoeff);
    this->OpWBound = SGPP::op_factory::createOperationHestonW(*this->BoundGrid, *this->wCoeff);
    this->OpWInner = SGPP::op_factory::createOperationHestonW(*this->InnerGrid, *this->wCoeff);
    this->OpZBound = SGPP::op_factory::createOperationHestonZ(*this->BoundGrid, *this->zCoeff);
    this->OpZInner = SGPP::op_factory::createOperationHestonZ(*this->InnerGrid, *this->zCoeff);
    this->OpGBound = SGPP::op_factory::createOperationHestonGLog(*this->BoundGrid, *this->gCoeff);
    this->OpGInner = SGPP::op_factory::createOperationHestonGLog(*this->InnerGrid, *this->gCoeff);
    this->OpDBound = SGPP::op_factory::createOperationHestonDLog(*this->BoundGrid, *this->dCoeff);
    this->OpDInner = SGPP::op_factory::createOperationHestonDLog(*this->InnerGrid, *this->dCoeff);
    this->OpFBound = SGPP::op_factory::createOperationHestonFLog(*this->BoundGrid, *this->fCoeff);
    this->OpFInner = SGPP::op_factory::createOperationHestonFLog(*this->InnerGrid, *this->fCoeff);
  } else {
    // Log-transformed case
    // Build coefficients
    buildBCoefficientsLogTransform();
    buildCCoefficientsLogTransform();
    buildDCoefficientsLogTransform();
    buildECoefficientsLogTransform();
    buildFCoefficientsLogTransform();
    buildGCoefficientsLogTransform();
    buildHCoefficientsLogTransform();
    buildKCoefficientsLogTransform();

    // Create operators
    this->OpBBound = SGPP::op_factory::createOperationHestonBLog(*this->BoundGrid, *this->bCoeff);
    this->OpBInner = SGPP::op_factory::createOperationHestonBLog(*this->InnerGrid, *this->bCoeff);
    this->OpCBound = SGPP::op_factory::createOperationHestonCLog(*this->BoundGrid, *this->cCoeff);
    this->OpCInner = SGPP::op_factory::createOperationHestonCLog(*this->InnerGrid, *this->cCoeff);
    this->OpDBound = SGPP::op_factory::createOperationHestonDLog(*this->BoundGrid, *this->dCoeff);
    this->OpDInner = SGPP::op_factory::createOperationHestonDLog(*this->InnerGrid, *this->dCoeff);
    this->OpEBound = SGPP::op_factory::createOperationHestonELog(*this->BoundGrid, *this->eCoeff);
    this->OpEInner = SGPP::op_factory::createOperationHestonELog(*this->InnerGrid, *this->eCoeff);
    this->OpFBound = SGPP::op_factory::createOperationHestonFLog(*this->BoundGrid, *this->fCoeff);
    this->OpFInner = SGPP::op_factory::createOperationHestonFLog(*this->InnerGrid, *this->fCoeff);
    this->OpGBound = SGPP::op_factory::createOperationHestonGLog(*this->BoundGrid, *this->gCoeff);
    this->OpGInner = SGPP::op_factory::createOperationHestonGLog(*this->InnerGrid, *this->gCoeff);
    this->OpHBound = SGPP::op_factory::createOperationHestonHLog(*this->BoundGrid, *this->hCoeff);
    this->OpHInner = SGPP::op_factory::createOperationHestonHLog(*this->InnerGrid, *this->hCoeff);
    this->OpKBound = SGPP::op_factory::createOperationHestonKLog(*this->BoundGrid, &(this->kCoeff));
    this->OpKInner = SGPP::op_factory::createOperationHestonKLog(*this->InnerGrid, &(this->kCoeff));
  }

  // Create operations that are needed by both the log-transformed and Cartesian case
  this->OpAInner = SGPP::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);
  this->OpABound = SGPP::op_factory::createOperationLTwoDotProduct(*this->BoundGrid);

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

HestonParabolicPDESolverSystemEuroAmer::~HestonParabolicPDESolverSystemEuroAmer() {
  delete this->dCoeff;
  delete this->eCoeff;
  delete this->fCoeff;
  delete this->gCoeff;
  delete this->zCoeff;
  delete this->bCoeff;
  delete this->cCoeff;
  delete this->hCoeff;
  delete this->xCoeff;
  delete this->yCoeff;
  delete this->wCoeff;
  delete4dEqualDimSizeArray(this->HestonAlgoDims.size(), &(this->kCoeff));

  if (!b_log_transform) {
    delete this->OpXBound;
    delete this->OpXInner;
    delete this->OpYBound;
    delete this->OpYInner;
    delete this->OpWBound;
    delete this->OpWInner;
    delete this->OpZBound;
    delete this->OpZInner;
    delete this->OpGBound;
    delete this->OpGInner;
    delete this->OpDBound;
    delete this->OpDInner;
    delete this->OpFBound;
    delete this->OpFInner;
  } else {
    delete this->OpBBound;
    delete this->OpBInner;
    delete this->OpCBound;
    delete this->OpCInner;
    delete this->OpDBound;
    delete this->OpDInner;
    delete this->OpEBound;
    delete this->OpEInner;
    delete this->OpFBound;
    delete this->OpFInner;
    delete this->OpGBound;
    delete this->OpGInner;
    delete this->OpHBound;
    delete this->OpHInner;
    delete this->OpKBound;
    delete this->OpKInner;
  }

  delete this->OpAInner;
  delete this->OpABound;

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
  delete this->oldGridStorage;
  delete this->secondGridStorage;
}

void HestonParabolicPDESolverSystemEuroAmer::applyLOperatorComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the riskfree rate
  if (this->r != 0.0) {
    this->OpABound->mult(alpha, temp);
    result.axpy((-1.0) * this->r, temp);
  }

  if (this->b_log_transform) {
    // Log-transformed coordinates

    // Apply the B operator
    this->OpBBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the C operator
    this->OpCBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the D operator
    this->OpDBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the E operator
    this->OpEBound->mult(alpha, temp);
    result.add(temp);

    // Apply the F operator
    this->OpFBound->mult(alpha, temp);
    result.add(temp);

    // Apply the G operator
    this->OpGBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the H operator
    this->OpHBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the K operator
    this->OpKBound->mult(alpha, temp);
    result.sub(temp);
  } else {
    // Cartesian coordinates
    // Only a single asset is supported.

    // Apply the X operator
    this->OpXBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the Y operator
    this->OpYBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the G operator
    this->OpGBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the D operator
    this->OpDBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the F operator
    this->OpFBound->mult(alpha, temp);
    result.add(temp);

    // Apply the W operator
    this->OpWBound->mult(alpha, temp);
    result.sub(temp);

    // Apply the Z operator
    this->OpZBound->mult(alpha, temp);
    result.add(temp);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::applyLOperatorInner(SGPP::base::DataVector& alpha,
                                                                 SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the riskfree rate
  if (this->r != 0.0) {
    this->OpAInner->mult(alpha, temp);
    result.axpy((-1.0) * this->r, temp);
  }

  if (this->b_log_transform) {
    // Log-transformed coordinates

    // Apply the B method
    this->OpBInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the C method
    this->OpCInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the D method
    this->OpDInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the E method
    this->OpEInner->mult(alpha, temp);
    result.add(temp);

    // Apply the F method
    this->OpFInner->mult(alpha, temp);
    result.add(temp);

    // Apply the G method
    this->OpGInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the H method
    this->OpHInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the K method
    this->OpKInner->mult(alpha, temp);
    result.sub(temp);
  } else {
    // Cartesian coordinates

    // Apply the X method
    this->OpXInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the Y method
    this->OpYInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the G method
    this->OpGInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the D method
    this->OpDInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the F method
    this->OpFInner->mult(alpha, temp);
    result.add(temp);

    // Apply the W method
    this->OpWInner->mult(alpha, temp);
    result.sub(temp);

    // Apply the Z method
    this->OpZInner->mult(alpha, temp);
    result.add(temp);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::applyMassMatrixComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpABound->mult(alpha, temp);

  result.add(temp);
}

void HestonParabolicPDESolverSystemEuroAmer::applyMassMatrixInner(SGPP::base::DataVector& alpha,
                                                                  SGPP::base::DataVector& result) {
  SGPP::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpAInner->mult(alpha, temp);

  result.add(temp);
}

void HestonParabolicPDESolverSystemEuroAmer::finishTimestep() {
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
}

void HestonParabolicPDESolverSystemEuroAmer::coarsenAndRefine(bool isLastTimestep) {
  // add number of Gridpoints
  this->numSumGridpointsInner += this->InnerGrid->getSize();
  this->numSumGridpointsComplete += this->BoundGrid->getSize();

  if (this->useCoarsen == true && isLastTimestep == false) {
    ///////////////////////////////////////////////////
    // Start integrated refinement & coarsening
    ///////////////////////////////////////////////////

    size_t originalGridSize = this->BoundGrid->getStorage().size();

    // Coarsen the grid
    std::unique_ptr<SGPP::base::GridGenerator> myGenerator = this->BoundGrid->createGridGenerator();

    if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numRefines = myGenerator->getNumberOfRefinablePoints();
      SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(
          this->alpha_complete, numRefines, this->refineThreshold);

      if (this->refineMode == "maxLevel") {
        myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage().size());
      }

      if (this->refineMode == "classic") {
        myGenerator->refine(myRefineFunc);
        this->alpha_complete->resizeZero(this->BoundGrid->getStorage().size());
      }

      delete myRefineFunc;
    }

    if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
      size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
      SGPP::base::SurplusCoarseningFunctor* myCoarsenFunctor =
          new SGPP::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen,
                                                   this->coarsenThreshold);
      myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
      delete myCoarsenFunctor;
    }

    ///////////////////////////////////////////////////
    // End integrated refinement & coarsening
    ///////////////////////////////////////////////////

    // rebuild the inner grid + coefficients
    this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                                   &this->InnerGrid, &this->alpha_inner);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::startTimestep() {
#ifndef NOBOUNDARYDISCOUNT

  // Adjust the boundaries with the riskfree rate
  if (this->r != 0.0) {
    if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
      SGPP::base::GridStorage* storage = &this->BoundGrid->getStorage();

      for (size_t i = 0; i < storage->size(); i++) {
        SGPP::base::GridIndex* curPoint = (*storage)[i];

        SGPP::base::DataVector pointCoords(storage->dim());
        curPoint->getCoords(pointCoords);

        // Discount the boundary points
        if (!curPoint->isInnerPoint()) {
          this->alpha_complete->set(
              i, this->alpha_complete->get(i) * exp(((-1.0) * (this->r * this->TimestepSize))));
        }
      }
    }
  }

#endif
}

void HestonParabolicPDESolverSystemEuroAmer::buildGCoefficients() {
  float_t volvol = volvols->get(0);
  float_t rho = this->hMatrix->get(0, 1);
  float_t kappa = this->kappas->get(0);

  this->gCoeff->setAll(0.0);
  this->gCoeff->set(1, (kappa) + rho * volvol);
}

void HestonParabolicPDESolverSystemEuroAmer::buildFCoefficients() {
  float_t volvol = volvols->get(0);
  float_t kappa = this->kappas->get(0);
  float_t theta = this->thetas->get(0);

  this->fCoeff->setAll(0.0);
  this->fCoeff->set(1, kappa * theta - 0.5 * pow(volvol, 2.0));
}

void HestonParabolicPDESolverSystemEuroAmer::buildDCoefficients() {
  float_t volvol = volvols->get(0);

  this->dCoeff->setAll(0.0);
  this->dCoeff->set(1, (0.5) * pow(volvol, 2.0));
}

void HestonParabolicPDESolverSystemEuroAmer::buildXCoefficients() {
  this->xCoeff->setAll(0.0);
  this->xCoeff->set(0, 1, 1.0);
}

void HestonParabolicPDESolverSystemEuroAmer::buildWCoefficients() {
  float_t volvol = volvols->get(0);
  float_t rho = this->hMatrix->get(0, 1);

  this->wCoeff->setAll(0.0);
  this->wCoeff->set(0, 1, rho * volvol);
}

void HestonParabolicPDESolverSystemEuroAmer::buildZCoefficients() {
  this->zCoeff->setAll(0.0);
  this->zCoeff->set(0, this->r);
}

void HestonParabolicPDESolverSystemEuroAmer::buildYCoefficients() {
  this->yCoeff->setAll(0.0);
  this->yCoeff->set(0, 1, 0.5);
}

void HestonParabolicPDESolverSystemEuroAmer::buildBCoefficientsLogTransform() {
  this->bCoeff->setAll(0.0);

  for (size_t i = 0; i < this->nAssets; i++) {
    this->bCoeff->set(2 * i, 2 * i + 1, 0.5);  // Coefficient is a constant of 0.5
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildCCoefficientsLogTransform() {
  float_t volvol, rho;
  this->cCoeff->setAll(0.0);

  for (size_t i = 0; i < this->nAssets; i++) {
    // Determine and set the coefficient for this asset
    volvol = volvols->get(i);
    rho = this->hMatrix->get(2 * i, 2 * i + 1);
    this->cCoeff->set(2 * i, 2 * i + 1, volvol * rho);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildDCoefficientsLogTransform() {
  this->dCoeff->setAll(0.0);

  for (size_t i = 0; i < this->nAssets; i++) {
    // Determine and set the coefficient for this asset
    float_t volvol = volvols->get(i);
    this->dCoeff->set(2 * i + 1, (0.5) * pow(volvol, 2.0));
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildECoefficientsLogTransform() {
  this->eCoeff->setAll(0.0);

  for (size_t i = 0; i < this->nAssets; i++) {
    this->eCoeff->set(2 * i, this->r);  // Coefficient is independent of which asset it is.
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildFCoefficientsLogTransform() {
  float_t theta, kappa, volvol;
  this->fCoeff->setAll(0.0);

  for (size_t i = 0; i < nAssets; i++) {
    theta = this->thetas->get(i);
    kappa = this->kappas->get(i);
    volvol = this->volvols->get(i);
    this->fCoeff->set(2 * i + 1, kappa * theta - 0.5 * pow(volvol, 2.0));
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildGCoefficientsLogTransform() {
  this->gCoeff->setAll(0.0);
  float_t kappa;

  for (size_t i = 0; i < nAssets; i++) {
    kappa = this->kappas->get(i);
    this->gCoeff->set(2 * i + 1, kappa);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildHCoefficientsLogTransform() {
  this->hCoeff->setAll(0.0);

  for (size_t i = 0; i < nAssets; i++) {
    this->hCoeff->set(2 * i, 2 * i + 1, 0.5);
  }
}

void HestonParabolicPDESolverSystemEuroAmer::buildKCoefficientsLogTransform() {
  float_t rho;
  setAll4dEqualDimSizeArray(this->HestonAlgoDims.size(), &(this->kCoeff), 0.0);

  for (size_t i = 0; i < nAssets; i++) {
    for (size_t j = i + 1; j < nAssets; j++) {
      rho = this->hMatrix->get(2 * i, 2 * j);
      this->kCoeff[2 * i][2 * i + 1][2 * j][2 * j + 1] = rho;
    }
  }
}

void HestonParabolicPDESolverSystemEuroAmer::create4dEqualDimSizeArray(size_t dimSize,
                                                                       float_t***** array) {
  (*array) = reinterpret_cast<float_t****>(calloc(dimSize, sizeof(float_t***)));

  for (size_t i = 0; i < dimSize; i++) {
    (*array)[i] = reinterpret_cast<float_t***>(calloc(dimSize, sizeof(float_t**)));

    for (size_t j = 0; j < dimSize; j++) {
      (*array)[i][j] = reinterpret_cast<float_t**>(calloc(dimSize, sizeof(float_t*)));

      for (size_t k = 0; k < dimSize; k++) {
        (*array)[i][j][k] = reinterpret_cast<float_t*>(calloc(dimSize, sizeof(float_t)));

        for (size_t m = 0; m < dimSize; m++) {
          (*array)[i][j][k][m] = 0.0;
        }
      }
    }
  }
}

void HestonParabolicPDESolverSystemEuroAmer::delete4dEqualDimSizeArray(size_t dimSize,
                                                                       float_t***** array) {
  for (size_t i = 0; i < dimSize; i++) {
    for (size_t j = 0; j < dimSize; j++) {
      for (size_t k = 0; k < dimSize; k++) {
        delete[](*array)[i][j][k];
      }

      delete[](*array)[i][j];
    }

    delete[](*array)[i];
  }

  delete[](*array);
}

void HestonParabolicPDESolverSystemEuroAmer::setAll4dEqualDimSizeArray(size_t dimSize,
                                                                       float_t***** array,
                                                                       float_t value) {
  for (size_t i = 0; i < dimSize; i++) {
    for (size_t j = 0; j < dimSize; j++) {
      for (size_t k = 0; k < dimSize; k++) {
        for (size_t m = 0; m < dimSize; m++) {
          (*array)[i][j][k][m] = value;
        }
      }
    }
  }
}
}  // namespace finance
}  // namespace SGPP
