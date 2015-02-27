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

#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    ModifiedBlackScholesParabolicPDESolverSystem::ModifiedBlackScholesParabolicPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& mu,
        SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho, float_t r, float_t TimestepSize, std::string OperationMode,
        bool bLogTransform, bool useCoarsen, float_t coarsenThreshold, std::string adaptSolveMode,
        int numCoarsenPoints, float_t refineThreshold, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int dim_HW)
      : BlackScholesParabolicPDESolverSystem(SparseGrid,
                                             alpha,
                                             mu,
                                             sigma,
                                             rho,
                                             r,
                                             TimestepSize,
                                             OperationMode,
                                             0.0,
                                             "nothing",
                                             bLogTransform,
                                             useCoarsen,
                                             coarsenThreshold,
                                             adaptSolveMode,
                                             numCoarsenPoints,
                                             refineThreshold,
                                             refineMode,
                                             refineMaxLevel) {
      this->OpFBound = SGPP::op_factory::createOperationLF(*this->BoundGrid);
      this->dim_r = dim_HW;
      this->variableDiscountFactor = new VariableDiscountFactor(SparseGrid.getStorage(), dim_HW);
    }

    void ModifiedBlackScholesParabolicPDESolverSystem::multiplyrBSHW(SGPP::base::DataVector& updateVector) {
      float_t tmp;

      for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++) {
        //std::string coords = (*storage)[i]->getCoordsStringBB(*this->myBoundingBox);
        std::string coords = this->BoundGrid->getStorage()->get(i)->getCoordsStringBB(*this->BoundGrid->getBoundingBox());
        std::stringstream coordsStream(coords);
        float_t dblFuncValues[2];

        for (size_t j = 0; j < 2; j++) {
          coordsStream >> tmp;
          dblFuncValues[j] = tmp;
        }

        // std::cout<< dblFuncValues[1]<< std::endl;
        //updateVector.set(i, updateVector.get(i)* dblFuncValues[1]);
        updateVector.set(i, updateVector.get(i)* dblFuncValues[this->dim_r]);
      }
    }

    ModifiedBlackScholesParabolicPDESolverSystem::~ModifiedBlackScholesParabolicPDESolverSystem() {
      delete this->OpFBound;
      delete this->variableDiscountFactor;
    }

    void ModifiedBlackScholesParabolicPDESolverSystem::applyLOperator(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the riskfree rate
      if (this->r != 0.0) {
        this->OpLTwoBound->mult(alpha, temp);
        result.axpy((-1.0)*this->r, temp);
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
      SGPP::base::DataVector factor(this->alpha_complete->getSize());
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
      if (this->useCoarsen == true) { //  && isLastTimestep == false)  // do it always as mostly only 1 timestep is executed
        ///////////////////////////////////////////////////
        // Start integrated refinement & coarsening
        ///////////////////////////////////////////////////

        size_t originalGridSize = this->BoundGrid->getStorage()->size();

        // Coarsen the grid
        SGPP::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

        //std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
        //std::cout << "Grid Size: " << originalGridSize << std::endl;

        if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
          size_t numRefines = myGenerator->getNumberOfRefinablePoints();
          SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);

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

        if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
          size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
          SGPP::base::SurplusCoarseningFunctor* myCoarsenFunctor = new SGPP::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
          myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
          delete myCoarsenFunctor;
        }

        delete myGenerator;

        ///////////////////////////////////////////////////
        // End integrated refinement & coarsening
        ///////////////////////////////////////////////////
      }
    }

    void ModifiedBlackScholesParabolicPDESolverSystem::startTimestep() {
      SGPP::base::DataVector factor(this->alpha_complete->getSize());
      // Adjust the boundaries with the riskfree rate

      this->variableDiscountFactor->getDiscountFactor(factor, this->TimestepSize);

      if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
        this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete, factor);
      }

    }
  }
}