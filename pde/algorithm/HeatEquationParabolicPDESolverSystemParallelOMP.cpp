/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp"
#include "base/exception/algorithm_exception.hpp"

#include "pde/algorithm/StdUpDown.hpp"
#include "pde/algorithm/UpDownOneOpDim.hpp"

#include "pde/operation/PdeOpFactory.hpp"
using namespace sg::op_factory;

#ifdef _OPENMP
#include "omp.h"
#endif

namespace sg {
  namespace pde {

    HeatEquationParabolicPDESolverSystemParallelOMP::HeatEquationParabolicPDESolverSystemParallelOMP(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode) {
      this->a = a;
      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;
      this->InnerGrid = NULL;
      this->alpha_inner = NULL;

      this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new sg::base::DirichletGridConverter();

      this->OpLaplaceBound = createOperationLaplace(SparseGrid);
      this->OpMassBound = sg::op_factory::createOperationLTwoDotProduct(SparseGrid);

      // create the inner grid
      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

      //Create needed operations, on inner grid
      this->OpLaplaceInner = createOperationLaplace(*this->InnerGrid);
      this->OpMassInner = sg::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);

      // right hand side if System
      this->rhs = NULL;
    }

    HeatEquationParabolicPDESolverSystemParallelOMP::~HeatEquationParabolicPDESolverSystemParallelOMP() {
      delete this->OpLaplaceBound;
      delete this->OpMassBound;
      delete this->OpLaplaceInner;
      delete this->OpMassInner;

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
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      sg::base::DataVector temp(alpha.getSize());

      ((StdUpDown*)(this->OpMassBound))->multParallelBuildingBlock(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      sg::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);

      std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
      size_t nDims = algoDims.size();
#ifdef _OPENMP
      omp_lock_t Mutex;
      omp_init_lock(&Mutex);
#endif

      // Apply Laplace, parallel in Dimensions
      for (size_t i = 0; i < nDims; i++) {
        #pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
        {
          sg::base::DataVector myResult(result.getSize());

          /// @todo (heinecke) discuss methods in order to avoid this cast
          ((UpDownOneOpDim*)(this->OpLaplaceBound))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

          // semaphore
#ifdef _OPENMP
          omp_set_lock(&Mutex);
#endif
          temp.add(myResult);
#ifdef _OPENMP
          omp_unset_lock(&Mutex);
#endif
        }
      }

      #pragma omp taskwait

#ifdef _OPENMP
      omp_destroy_lock(&Mutex);
#endif

      result.axpy((-1.0)*this->a, temp);
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      sg::base::DataVector temp(alpha.getSize());

      ((StdUpDown*)(this->OpMassInner))->multParallelBuildingBlock(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      sg::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);

      std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
      size_t nDims = algoDims.size();
#ifdef _OPENMP
      omp_lock_t Mutex;
      omp_init_lock(&Mutex);
#endif

      // Apply Laplace, parallel in Dimensions
      for (size_t i = 0; i < nDims; i++) {
        #pragma omp task firstprivate(i) shared(alpha, temp, result, algoDims)
        {
          sg::base::DataVector myResult(result.getSize());

          /// @todo (heinecke) discuss methods in order to avoid this cast
          ((UpDownOneOpDim*)(this->OpLaplaceInner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

          // semaphore
#ifdef _OPENMP
          omp_set_lock(&Mutex);
#endif
          temp.add(myResult);
#ifdef _OPENMP
          omp_unset_lock(&Mutex);
#endif
        }
      }

      #pragma omp taskwait

#ifdef _OPENMP
      omp_destroy_lock(&Mutex);
#endif

      result.axpy((-1.0)*this->a, temp);
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::finishTimestep() {
      // Replace the inner coefficients on the boundary grid
      this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::coarsenAndRefine(bool isLastTimestep) {
    }


    void HeatEquationParabolicPDESolverSystemParallelOMP::startTimestep() {
    }

    void HeatEquationParabolicPDESolverSystemParallelOMP::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      if (this->tOperationMode == "ExEul") {
        applyMassMatrixInner(alpha, result);
      } else if (this->tOperationMode == "ImEul") {
        sg::base::DataVector temp(result.getSize());
        sg::base::DataVector temp2(result.getSize());

        #pragma omp parallel shared(alpha, result, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (alpha, temp)
            {
              applyMassMatrixInner(alpha, temp);
            }

            #pragma omp task shared (alpha, temp2)
            {
              applyLOperatorInner(alpha, temp2);
            }

            #pragma omp taskwait
          }
        }

        result.add(temp);
        result.axpy((-1.0)*this->TimestepSize, temp2);
      } else if (this->tOperationMode == "CrNic") {
        sg::base::DataVector temp(result.getSize());
        sg::base::DataVector temp2(result.getSize());

        #pragma omp parallel shared(alpha, result, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (alpha, temp)
            {
              applyMassMatrixInner(alpha, temp);
            }

            #pragma omp task shared (alpha, temp2)
            {
              applyLOperatorInner(alpha, temp2);
            }

            #pragma omp taskwait
          }
        }

        result.add(temp);
        result.axpy((-0.5)*this->TimestepSize, temp2);
      } else {
        throw new sg::base::algorithm_exception(" HeatEquationParabolicPDESolverSystemParallelOMP::mult : An unknown operation mode was specified!");
      }
    }

    sg::base::DataVector* HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS() {
      sg::base::DataVector rhs_complete(this->alpha_complete->getSize());

      if (this->tOperationMode == "ExEul") {
        rhs_complete.setAll(0.0);

        sg::base::DataVector temp(rhs_complete.getSize());
        sg::base::DataVector temp2(rhs_complete.getSize());
        sg::base::DataVector myAlpha(*this->alpha_complete);

        #pragma omp parallel shared(myAlpha, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (myAlpha, temp)
            {
              applyMassMatrixComplete(myAlpha, temp);
            }

            #pragma omp task shared (myAlpha, temp2)
            {
              applyLOperatorComplete(myAlpha, temp2);
            }

            #pragma omp taskwait
          }
        }

        rhs_complete.add(temp);
        rhs_complete.axpy(this->TimestepSize, temp2);
      } else if (this->tOperationMode == "ImEul") {
        rhs_complete.setAll(0.0);

        applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
      } else if (this->tOperationMode == "CrNic") {
        rhs_complete.setAll(0.0);

        sg::base::DataVector temp(rhs_complete.getSize());
        sg::base::DataVector temp2(rhs_complete.getSize());
        sg::base::DataVector myAlpha(*this->alpha_complete);

        #pragma omp parallel shared(myAlpha, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (myAlpha, temp)
            {
              applyMassMatrixComplete(myAlpha, temp);
            }

            #pragma omp task shared (myAlpha, temp2)
            {
              applyLOperatorComplete(myAlpha, temp2);
            }

            #pragma omp taskwait
          }
        }

        rhs_complete.add(temp);
        rhs_complete.axpy((0.5)*this->TimestepSize, temp2);
      } else {
        throw new sg::base::algorithm_exception("HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode was specified!");
      }

      this->startTimestep();

      // Now apply the boundary ansatzfunctions to the inner ansatzfunctions
      sg::base::DataVector result_complete(this->alpha_complete->getSize());
      sg::base::DataVector alpha_bound(*this->alpha_complete);

      result_complete.setAll(0.0);

      this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

      // apply CG Matrix
      if (this->tOperationMode == "ExEul") {
        applyMassMatrixComplete(alpha_bound, result_complete);
      } else if (this->tOperationMode == "ImEul") {
        sg::base::DataVector temp(alpha_bound.getSize());
        sg::base::DataVector temp2(alpha_bound.getSize());

        #pragma omp parallel shared(alpha_bound, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (alpha_bound, temp)
            {
              applyMassMatrixComplete(alpha_bound, temp);
            }

            #pragma omp task shared (alpha_bound, temp2)
            {
              applyLOperatorComplete(alpha_bound, temp2);
            }

            #pragma omp taskwait
          }
        }

        result_complete.add(temp);
        result_complete.axpy((-1.0)*this->TimestepSize, temp2);
      } else if (this->tOperationMode == "CrNic") {
        sg::base::DataVector temp(alpha_bound.getSize());
        sg::base::DataVector temp2(alpha_bound.getSize());

        #pragma omp parallel shared(alpha_bound, temp, temp2)
        {
          #pragma omp single nowait
          {
            #pragma omp task shared (alpha_bound, temp)
            {
              applyMassMatrixComplete(alpha_bound, temp);
            }

            #pragma omp task shared (alpha_bound, temp2)
            {
              applyLOperatorComplete(alpha_bound, temp2);
            }

            #pragma omp taskwait
          }
        }

        result_complete.add(temp);
        result_complete.axpy((-0.5)*this->TimestepSize, temp2);
      } else {
        throw new sg::base::algorithm_exception("HeatEquationParabolicPDESolverSystemParallelOMP::generateRHS : An unknown operation mode was specified!");
      }

      rhs_complete.sub(result_complete);

      if (this->rhs != NULL) {
        delete this->rhs;
      }

      this->rhs = new sg::base::DataVector(this->alpha_inner->getSize());
      this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

      return this->rhs;
    }

  }
}
