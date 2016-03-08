// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>

#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/pde/algorithm/StdUpDown.hpp>
#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

using sgpp::op_factory::createOperationLaplace;

namespace sgpp {
namespace parallel {

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::
    PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(sgpp::base::Grid& SparseGrid,
                                                               sgpp::base::DataVector& rhs)
    : sgpp::pde::OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
  this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
  this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
}

PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::
    ~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI() {}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

#pragma omp parallel shared(alpha, result)
  {
#pragma omp single nowait
    {
      std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
      size_t nDims = algoDims.size();

      // Apply Laplace, parallel in Dimensions
      for (size_t i = 0; i < nDims; i++) {
        if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
          {
            sgpp::base::DataVector myResult(result.getSize());

            /// discuss methods in order to avoid this cast
            ((sgpp::pde::UpDownOneOpDim*)(this->Laplace_Inner.get()))
                ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

#pragma omp critical
            { result.add(myResult); }
          }
        }
      }

#pragma omp taskwait
    }
  }
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::applyLOperatorComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

#pragma omp parallel shared(alpha, result)
  {
#pragma omp single nowait
    {
      std::vector<size_t> algoDims = this->InnerGrid->getStorage().getAlgorithmicDimensions();
      size_t nDims = algoDims.size();

      // Apply Laplace, parallel in Dimensions
      for (size_t i = 0; i < nDims; i++) {
        if (i % myGlobalMPIComm->getNumRanks() == myGlobalMPIComm->getMyRank()) {
#pragma omp task firstprivate(i) shared(alpha, result, algoDims)
          {
            sgpp::base::DataVector myResult(result.getSize());

            /// discuss methods in order to avoid this cast
            ((sgpp::pde::UpDownOneOpDim*)(this->Laplace_Complete.get()))
                ->multParallelBuildingBlock(alpha, myResult, algoDims[i]);

#pragma omp critical
            { result.add(myResult); }
          }
        }
      }

#pragma omp taskwait
    }
  }
}

void PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::mult(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  // distribute the current grid coefficients
  // myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);

  this->applyLOperatorInner(alpha, result);

  // aggregate all results
  // myGlobalMPIComm->reduceGridCoefficientsOnRank0(result);
  myGlobalMPIComm->reduceGridCoefficients(result);
}

sgpp::base::DataVector* PoissonEquationEllipticPDESolverSystemDirichletParallelMPI::generateRHS() {
  if (this->InnerGrid != NULL) {
    sgpp::base::DataVector alpha_tmp_complete(*(this->rhs));
    sgpp::base::DataVector rhs_tmp_complete(*(this->rhs));

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
    throw sgpp::base::algorithm_exception(
        "OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
  }

  return this->rhs_inner;
}
}  // namespace parallel
}  // namespace sgpp
