// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystemDirichlet.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <sstream>

namespace sgpp {
namespace pde {

OperationEllipticPDESolverSystemDirichlet::OperationEllipticPDESolverSystemDirichlet(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& rhs)
    : OperationEllipticPDESolverSystem(SparseGrid, rhs) {
  this->BoundaryUpdate = new sgpp::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->GridConverter = new sgpp::base::DirichletGridConverter();

  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->rhs, &(this->InnerGrid),
                                               &(this->rhs_inner));

  this->numGridpointsInner = this->InnerGrid->getSize();

  this->alpha_inner = nullptr;
}

OperationEllipticPDESolverSystemDirichlet::~OperationEllipticPDESolverSystemDirichlet() {
  delete this->alpha_inner;
  delete this->rhs_inner;
  delete this->InnerGrid;
  delete this->BoundaryUpdate;
  delete this->GridConverter;
}

void OperationEllipticPDESolverSystemDirichlet::mult(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result) {
  applyLOperatorInner(alpha, result);
}

sgpp::base::DataVector* OperationEllipticPDESolverSystemDirichlet::generateRHS() {
  if (this->InnerGrid != nullptr) {
    sgpp::base::DataVector alpha_tmp_complete(*(this->rhs));
    sgpp::base::DataVector rhs_tmp_complete(*(this->rhs));

    this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
    applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

    this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
    this->rhs_inner->mult(-1.0);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
  }

  return this->rhs_inner;
}

sgpp::base::DataVector* OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG() {
  if (this->InnerGrid != nullptr) {
    if (this->alpha_inner != nullptr) {
      delete this->alpha_inner;
    }

    this->alpha_inner = new sgpp::base::DataVector(this->InnerGrid->getSize());
    this->alpha_inner->setAll(0.0);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG : No inner grid "
        "exists!");
  }

  return this->alpha_inner;
}

void OperationEllipticPDESolverSystemDirichlet::getSolutionBoundGrid(
    sgpp::base::DataVector& Solution, sgpp::base::DataVector& SolutionInner) {
  Solution = *(this->rhs);
  this->GridConverter->updateBoundaryCoefs(Solution, SolutionInner);
}
}  // namespace pde
}  // namespace sgpp
