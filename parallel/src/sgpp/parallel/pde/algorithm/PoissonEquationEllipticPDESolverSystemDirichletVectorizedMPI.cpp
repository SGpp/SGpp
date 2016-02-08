// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/operation/ParallelOpFactory.hpp>
#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(
  SGPP::base::Grid& SparseGrid,
  SGPP::base::DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(
      SparseGrid, rhs) {
  // Create operations
  char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

  if (! strcmp(alg_selector, "X86SIMD")) {
    this->Laplace_Inner = SGPP::op_factory::createOperationLaplaceVectorized(
                            *this->InnerGrid, SGPP::parallel::X86SIMD);
    this->Laplace_Complete = SGPP::op_factory::createOperationLaplaceVectorized(
                               *this->BoundGrid, SGPP::parallel::X86SIMD);
#ifdef USEOCL
  } else if (! strcmp(alg_selector, "OCL")) {
    this->Laplace_Inner = SGPP::op_factory::createOperationLaplaceVectorized(
                            *this->InnerGrid, SGPP::parallel::OpenCL);
    this->Laplace_Complete = SGPP::op_factory::createOperationLaplaceVectorized(
                               *this->BoundGrid, SGPP::parallel::OpenCL);
#endif
  } else {
    throw SGPP::base::algorithm_exception("PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization was selected!");
  }
}

PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::~PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI() {
  delete this->Laplace_Complete;
  delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorInner(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  Laplace_Inner->mult(alpha, result);
}

void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorComplete(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  Laplace_Complete->mult(alpha, result);
}

}
}