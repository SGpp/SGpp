/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/operation/ParallelOpFactory.hpp"
#include "parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI.hpp"

#include "base/exception/algorithm_exception.hpp"

#include <cstring>

namespace sg {
  namespace parallel {

    PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
      // Create operations
      char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

      if (! strcmp(alg_selector, "X86SIMD")) {
        this->Laplace_Inner = sg::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, sg::parallel::X86SIMD);
        this->Laplace_Complete = sg::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, sg::parallel::X86SIMD);
#ifdef USEOCL
      } else if (! strcmp(alg_selector, "OCL")) {
        this->Laplace_Inner = sg::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, sg::parallel::OpenCL);
        this->Laplace_Complete = sg::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, sg::parallel::OpenCL);
#endif
      } else {
        throw sg::base::algorithm_exception("PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization was selected!");
      }
    }

    PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::~PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI() {
      delete this->Laplace_Complete;
      delete this->Laplace_Inner;
    }

    void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      Laplace_Inner->mult(alpha, result);
    }

    void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      Laplace_Complete->mult(alpha, result);
    }

  }
}
