// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace pde {

HeatEquationParabolicPDESolverSystem::HeatEquationParabolicPDESolverSystem(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, double a, double TimestepSize,
    std::string OperationMode) {
  this->a = a;
  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;
  this->InnerGrid = nullptr;
  this->alpha_inner = nullptr;

  this->BoundaryUpdate = new sgpp::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->GridConverter = new sgpp::base::DirichletGridConverter();

  this->OpLaplaceBound = op_factory::createOperationLaplace(SparseGrid);
  this->OpMassBound = sgpp::op_factory::createOperationLTwoDotProduct(SparseGrid);

  // create the inner grid
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                               &this->InnerGrid, &this->alpha_inner);

  // Create needed operations, on inner grid
  this->OpLaplaceInner = op_factory::createOperationLaplace(*this->InnerGrid);
  this->OpMassInner = sgpp::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);

  // right hand side if System
  this->rhs = new sgpp::base::DataVector(1);
}

HeatEquationParabolicPDESolverSystem::~HeatEquationParabolicPDESolverSystem() {
  delete this->OpLaplaceBound;
  delete this->OpMassBound;
  delete this->OpLaplaceInner;
  delete this->OpMassInner;

  delete this->BoundaryUpdate;
  delete this->GridConverter;

  if (this->InnerGrid != nullptr) {
    delete this->InnerGrid;
  }

  if (this->alpha_inner != nullptr) {
    delete this->alpha_inner;
  }

  delete this->rhs;
}

void HeatEquationParabolicPDESolverSystem::applyMassMatrixComplete(sgpp::base::DataVector& alpha,
                                                                   sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  // Apply the mass matrix
  this->OpMassBound->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystem::applyLOperatorComplete(sgpp::base::DataVector& alpha,
                                                                  sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  // Apply the laplace Operator rate
  this->OpLaplaceBound->mult(alpha, temp);
  result.axpy((-1.0) * this->a, temp);
}

void HeatEquationParabolicPDESolverSystem::applyMassMatrixInner(sgpp::base::DataVector& alpha,
                                                                sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  // Apply the mass matrix
  this->OpMassInner->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystem::applyLOperatorInner(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result) {
  result.setAll(0.0);

  sgpp::base::DataVector temp(alpha.getSize());

  // Apply the laplace Operator rate
  this->OpLaplaceInner->mult(alpha, temp);
  result.axpy((-1.0) * this->a, temp);
}

void HeatEquationParabolicPDESolverSystem::finishTimestep() {
  // Replace the inner coefficients on the boundary grid
  this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {}

void HeatEquationParabolicPDESolverSystem::startTimestep() {}
}  // namespace pde
}  // namespace sgpp
