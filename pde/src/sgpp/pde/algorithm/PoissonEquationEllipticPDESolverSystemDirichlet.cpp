// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

PoissonEquationEllipticPDESolverSystemDirichlet::PoissonEquationEllipticPDESolverSystemDirichlet(
    SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs)
    : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs) {
  this->Laplace_Complete = op_factory::createOperationLaplace(*this->BoundGrid).release();
  this->Laplace_Inner = op_factory::createOperationLaplace(*this->InnerGrid).release();
}

PoissonEquationEllipticPDESolverSystemDirichlet::
    ~PoissonEquationEllipticPDESolverSystemDirichlet() {
  delete this->Laplace_Complete;
  delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorInner(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  Laplace_Inner->mult(alpha, result);
}

void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorComplete(
    SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  Laplace_Complete->mult(alpha, result);
}
}  // namespace pde
}  // namespace SGPP
