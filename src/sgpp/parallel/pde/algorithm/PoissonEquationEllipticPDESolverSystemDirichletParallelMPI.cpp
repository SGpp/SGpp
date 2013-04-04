/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/tools/MPI/SGppMPITools.hpp"

#include "parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "pde/operation/PdeOpFactory.hpp"

#include "pde/algorithm/StdUpDown.hpp"
#include "pde/algorithm/UpDownOneOpDim.hpp"

using namespace sg::op_factory;

namespace sg {
  namespace parallel {

    PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : sg::pde::OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
      this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
      this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
    }

    PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI() {
      delete this->Laplace_Complete;
      delete this->Laplace_Inner;
    }

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      #pragma omp parallel shared(alpha, result)
      {
        #pragma omp single nowait
        {
          std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
          size_t nDims = algoDims.size();

          // Apply Laplace, parallel in Dimensions
          for (size_t i = 0; i < nDims; i++) {
            if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
              #pragma omp task firstprivate(i) shared(alpha, result, algoDims)
              {
                sg::base::DataVector myResult(result.getSize());

                /// @todo (heinecke) discuss methods in order to avoid this cast
                ((sg::pde::UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

                #pragma omp critical
                {
                  result.add(myResult);
                }
              }
            }
          }

          #pragma omp taskwait
        }
      }

    }

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      #pragma omp parallel shared(alpha, result)
      {
        #pragma omp single nowait
        {
          std::vector<size_t> algoDims = this->InnerGrid->getStorage()->getAlgorithmicDimensions();
          size_t nDims = algoDims.size();

          // Apply Laplace, parallel in Dimensions
          for (size_t i = 0; i < nDims; i++) {
            if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
              #pragma omp task firstprivate(i) shared(alpha, result, algoDims)
              {

                sg::base::DataVector myResult(result.getSize());

                /// @todo (heinecke) discuss methods in order to avoid this cast
                ((sg::pde::UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

                #pragma omp critical
                {
                  result.add(myResult);
                }
              }
            }
          }

          #pragma omp taskwait
        }
      }
    }

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      // distribute the current grid coefficients
      //myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

      this->applyLOperatorInner(alpha, result);

      // aggregate all results
      //myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
      myGlobalMPIComm->reduceGridCoefficients(result);
    }

    sg::base::DataVector* PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::generateRHS() {
      if (this->InnerGrid != NULL) {
        sg::base::DataVector alpha_tmp_complete(*(this->rhs));
        sg::base::DataVector rhs_tmp_complete(*(this->rhs));

        this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
        // distribute the current grid coefficients
        myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha_tmp_complete);

        applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

        // aggregate all results
        myGlobalMPIComm->reduceGridCoefficientsOnRank0(rhs_tmp_complete);

        this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
        this->rhs_inner->mult(-1.0);
      } else {
        myGlobalMPIComm->Abort();
        throw new sg::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
      }

      return this->rhs_inner;
    }

  }
}
