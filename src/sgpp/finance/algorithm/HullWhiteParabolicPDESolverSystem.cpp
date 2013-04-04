/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi (qic@in.tum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "finance/algorithm/HullWhiteParabolicPDESolverSystem.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/Grid.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "finance/operation/FinanceOpFactory.hpp"

#include <cmath>

namespace sg {
  namespace finance {

    HullWhiteParabolicPDESolverSystem::HullWhiteParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double sigma,
        double theta, double a, double TimestepSize, std::string OperationMode,
        bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
        int numCoarsenPoints, double refineThreshold, std::string refineMode, sg::base::GridIndex::level_type refineMaxLevel,
        int dim_HW) {
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;

      this->alpha_complete_old = new sg::base::DataVector(this->alpha_complete->getSize());
      this->alpha_complete_tmp = new sg::base::DataVector(this->alpha_complete->getSize());
      this->oldGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());
      this->secondGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());

      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->TimestepSize_old = TimestepSize;
      this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->variableDiscountFactor = new VariableDiscountFactor(SparseGrid.getStorage(), dim_HW);
      this->a = a;
      this->theta = theta;
      this->sigma = sigma;
      this->HWalgoDims = this->BoundGrid->getAlgorithmicDimensions();

      // Create needed operations, on boundary grid
      this->OpBBound = sg::op_factory::createOperationLB(*this->BoundGrid);
      this->OpDBound = sg::op_factory::createOperationLD(*this->BoundGrid);
      this->OpEBound = sg::op_factory::createOperationLE(*this->BoundGrid);
      this->OpFBound = sg::op_factory::createOperationLF(*this->BoundGrid);

      // Create operations, independent bLogTransform
      this->OpLTwoBound = sg::op_factory::createOperationLTwoDotProduct(*this->BoundGrid);

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
      this->dim_r = dim_HW;
    }

    HullWhiteParabolicPDESolverSystem::~HullWhiteParabolicPDESolverSystem() {
      delete this->OpBBound;
      delete this->OpDBound;
      delete this->OpEBound;
      delete this->OpFBound;
      delete this->OpLTwoBound;
      delete this->BoundaryUpdate;
      delete this->variableDiscountFactor;

      if (this->rhs != NULL) {
        delete this->rhs;
      }

      delete this->alpha_complete_old;
      delete this->alpha_complete_tmp;
      delete this->oldGridStorage;
      delete this->secondGridStorage;
    }

    void HullWhiteParabolicPDESolverSystem::applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      if (this->theta != 0.0) {
        this->OpBBound->mult(alpha, temp);
        result.axpy(1.0 * this->theta, temp);
      }

      if (this->sigma != 0.0) {
        this->OpEBound->mult(alpha, temp);
        result.axpy((-1.0 / 2.0)*pow((this->sigma), 2.0), temp);
      }

      if (this->a != 0.0) {
        this->OpFBound->mult(alpha, temp);
        result.axpy((-1.0)*this->a, temp);
      }


      this->OpDBound->mult(alpha, temp);
      result.sub(temp);
    }

    void HullWhiteParabolicPDESolverSystem::applyMassMatrix(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoBound->mult(alpha, temp);

      result.add(temp);
    }

    void HullWhiteParabolicPDESolverSystem::finishTimestep() {
      sg::base::DataVector factor(this->alpha_complete->getSize());
      // Adjust the boundaries with the riskfree rate
      this->variableDiscountFactor->getDiscountFactor(factor, this->TimestepSize);

      if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas") {
        this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete, factor);
      }

      // add number of Gridpoints
      this->numSumGridpointsInner += 0;
      this->numSumGridpointsComplete += this->BoundGrid->getSize();

    }

    void HullWhiteParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {
      if (this->useCoarsen == true) { // && isLastTimestep == false) // do it always as mostly only 1 timestep is executed
        ///////////////////////////////////////////////////
        // Start integrated refinement & coarsening
        ///////////////////////////////////////////////////

        size_t originalGridSize = this->BoundGrid->getStorage()->size();

        // Coarsen the grid
        sg::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

        //std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
        //std::cout << "Grid Size: " << originalGridSize << std::endl;

        if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
          size_t numRefines = myGenerator->getNumberOfRefinablePoints();
          sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);

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
          sg::base::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
          myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
          delete myCoarsenFunctor;
        }

        delete myGenerator;

        ///////////////////////////////////////////////////
        // End integrated refinement & coarsening
        ///////////////////////////////////////////////////
      }
    }

    void HullWhiteParabolicPDESolverSystem::startTimestep() {
      sg::base::DataVector factor(this->alpha_complete->getSize());
      // Adjust the boundaries with the riskfree rate
      this->variableDiscountFactor->getDiscountFactor(factor, this->TimestepSize);

      if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
        this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete, factor);
      }

    }

  }
}
