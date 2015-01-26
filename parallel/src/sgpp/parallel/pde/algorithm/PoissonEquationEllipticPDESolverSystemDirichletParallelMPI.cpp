/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>

#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/pde/algorithm/StdUpDown.hpp>
#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

using namespace SGPP::op_factory;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs) : SGPP::pde::OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
      this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
      this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
    }

    PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI() {
      delete this->Laplace_Complete;
      delete this->Laplace_Inner;
    }

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
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
                SGPP::base::DataVector myResult(result.getSize());

                /// @todo (heinecke) discuss methods in order to avoid this cast
                ((SGPP::pde::UpDownOneOpDim*)(this->Laplace_Inner))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

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

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
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

                SGPP::base::DataVector myResult(result.getSize());

                /// @todo (heinecke) discuss methods in order to avoid this cast
                ((SGPP::pde::UpDownOneOpDim*)(this->Laplace_Complete))->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

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

    void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      // distribute the current grid coefficients
      //myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

      this->applyLOperatorInner(alpha, result);

      // aggregate all results
      //myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
      myGlobalMPIComm->reduceGridCoefficients(result);
    }

    SGPP::base::DataVector* PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::generateRHS() {
      if (this->InnerGrid != NULL) {
        SGPP::base::DataVector alpha_tmp_complete(*(this->rhs));
        SGPP::base::DataVector rhs_tmp_complete(*(this->rhs));

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
        throw new SGPP::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
      }

      return this->rhs_inner;
    }

  }
}
