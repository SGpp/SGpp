// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/ModifiedBlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/finance/tools/VariableDiscountFactor.hpp>
#include <sgpp/finance/operation/FinanceOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <string>

namespace sgpp {
namespace finance {

ModifiedBlackScholesParabolicPDESolverSystem::ModifiedBlackScholesParabolicPDESolverSystem(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, sgpp::base::DataVector& mu,
    sgpp::base::DataVector& sigma, sgpp::base::DataMatrix& rho, double r, double TimestepSize,
    std::string OperationMode, bool bLogTransform, bool useCoarsen, double coarsenThreshold,
    std::string adaptSolveMode, int numCoarsenPoints, double refineThreshold,
    std::string refineMode, sgpp::base::GridPoint::level_type refineMaxLevel, int dim_HW)
    : BlackScholesParabolicPDESolverSystem(SparseGrid, alpha, mu, sigma, rho, r, TimestepSize,
                                           OperationMode, 0.0, "nothing", bLogTransform, useCoarsen,
                                           coarsenThreshold, adaptSolveMode, numCoarsenPoints,
                                           refineThreshold, refineMode, refineMaxLevel) {
  this->OpFBound = sgpp::op_factory::createOperationLF(*this->BoundGrid).release();
  this->dim_r = dim_HW;
  this->variableDiscountFactor = new VariableDiscountFactor(&SparseGrid.getStorage(), dim_HW);
}

void ModifiedBlackScholesParabolicPDESolverSystem::multiplyrBSHW(
    sgpp::base::DataVector& updateVector) {
  double tmp;

  for (size_t i = 0; i < this->BoundGrid->getSize(); i++) {
    // std::string coords = (*storage)[i].getCoordsStringBB(*this->myBoundingBox);
    std::string coords = this->BoundGrid->getStorage().getGridPoint(i).getCoordsStringBB(
        this->BoundGrid->getBoundingBox());
    std::stringstream coordsStream(coords);
    double dblFuncValues[2];

    for (size_t j = 0; j < 2; j++) {
      coordsStream >> tmp;
      dblFuncValues[j] = tmp;
    }

    // std::cout<< dblFuncValues[1]<< std::endl;
    // updateVector.set(i, updateVector.get(i)* dblFuncValues[1]);
    updateVector.set(i, updateVector.get(i) * dblFuncValues[this->dim_r]);
  }
}

ModifiedBlackScholesParabolicPDESolverSystem::~ModifiedBlackScholesParabolicPDESolverSystem() {
  delete this->OpFBound;
  delete this->variableDiscountFactor;
}

void ModifiedBlackScholesParabolicPDESolverSystem::applyLOperator(sgpp::base::DataVector& alpha,
                                                                  sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

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

  this->OpFBound->mult(alpha, temp);
  this->multiplyrBSHW(temp);
  result.add(temp);
}

void ModifiedBlackScholesParabolicPDESolverSystem::finishTimestep() {
  sgpp::base::DataVector factor(this->alpha_complete->getSize());
  // Adjust the boundaries with the riskfree rate
  this->variableDiscountFactor->getDiscountFactor(factor, this->TimestepSize);

  if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas") {
    this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete, factor);
  }

  // add number of Gridpoints
  this->numSumGridpointsInner += 0;
  this->numSumGridpointsComplete += this->BoundGrid->getSize();
}

void ModifiedBlackScholesParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {
  if (this->useCoarsen == true) {  //  && isLastTimestep == false)  // do it always as mostly only 1
                                   //  timestep is executed
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

void ModifiedBlackScholesParabolicPDESolverSystem::startTimestep() {
  sgpp::base::DataVector factor(this->alpha_complete->getSize());
  // Adjust the boundaries with the riskfree rate

  this->variableDiscountFactor->getDiscountFactor(factor, this->TimestepSize);

  if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
    this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete, factor);
  }
}
}  // namespace finance
}  // namespace sgpp
