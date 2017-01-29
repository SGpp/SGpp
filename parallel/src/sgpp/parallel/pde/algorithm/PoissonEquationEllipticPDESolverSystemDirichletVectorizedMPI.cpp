// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/operation/ParallelOpFactory.hpp>
#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {
namespace parallel {

PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::
    PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(sgpp::base::Grid& SparseGrid,
                                                                 sgpp::base::DataVector& rhs)
    : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
  // Create operations
  char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

  if (!strcmp(alg_selector, "X86SIMD")) {
    this->Laplace_Inner.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::X86SIMD));
    this->Laplace_Complete.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::X86SIMD));
#ifdef USEOCL
  } else if (!strcmp(alg_selector, "OCL")) {
    this->Laplace_Inner.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::OpenCL));
    this->Laplace_Complete.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::OpenCL));
#endif
  } else {
    throw sgpp::base::algorithm_exception(
        "PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::"
        "PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization "
        "was selected!");
  }
}

PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::
    ~PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI() {}

void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  Laplace_Inner->mult(alpha, result);
}

void PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::applyLOperatorComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  Laplace_Complete->mult(alpha, result);
}
}  // namespace parallel
}  // namespace sgpp
