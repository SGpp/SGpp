// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystemFreeBoundaries.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationEllipticPDESolverSystemFreeBoundaries::OperationEllipticPDESolverSystemFreeBoundaries(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& rhs)
    : OperationEllipticPDESolverSystem(SparseGrid, rhs) {}

OperationEllipticPDESolverSystemFreeBoundaries::~OperationEllipticPDESolverSystemFreeBoundaries() {}

void OperationEllipticPDESolverSystemFreeBoundaries::mult(sgpp::base::DataVector& alpha,
                                                          sgpp::base::DataVector& result) {
  applyLOperator(alpha, result);
}

sgpp::base::DataVector* OperationEllipticPDESolverSystemFreeBoundaries::generateRHS() {
  return this->rhs;
}
}  // namespace pde
}  // namespace sgpp
