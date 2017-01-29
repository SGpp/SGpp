// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemVectorizedMPI.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>
#include <string>

namespace sgpp {
namespace parallel {

HeatEquationParabolicPDESolverSystemVectorizedMPI::
    HeatEquationParabolicPDESolverSystemVectorizedMPI(sgpp::base::Grid& SparseGrid,
                                                      sgpp::base::DataVector& alpha, double a,
                                                      double TimestepSize,
                                                      std::string OperationMode) {
  this->a = a;
  this->tOperationMode = OperationMode;
  this->TimestepSize = TimestepSize;
  this->BoundGrid = &SparseGrid;
  this->alpha_complete = &alpha;
  this->InnerGrid = NULL;
  this->alpha_inner = NULL;

  this->BoundaryUpdate = new sgpp::base::DirichletUpdateVector(SparseGrid.getStorage());
  this->GridConverter = new sgpp::base::DirichletGridConverter();

  // create the inner grid
  this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete,
                                               &this->InnerGrid, &this->alpha_inner);

  // Create operations
  char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

  if (!strcmp(alg_selector, "X86SIMD")) {
    this->OpLaplaceInner.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::X86SIMD));
    this->OpLaplaceBound.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::X86SIMD));
    this->OpLTwoInner.reset(sgpp::op_factory::createOperationLTwoDotProductVectorized(
        *this->InnerGrid, sgpp::parallel::X86SIMD));
    this->OpLTwoBound.reset(sgpp::op_factory::createOperationLTwoDotProductVectorized(
        *this->BoundGrid, sgpp::parallel::X86SIMD));
    this->OpLTwoDotLaplaceInner.reset(sgpp::op_factory::createOperationLTwoDotLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::X86SIMD));
    this->OpLTwoDotLaplaceBound.reset(sgpp::op_factory::createOperationLTwoDotLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::X86SIMD));
#ifdef USEOCL
  } else if (!strcmp(alg_selector, "OCL")) {
    this->OpLaplaceInner.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::OpenCL));
    this->OpLaplaceBound.reset(sgpp::op_factory::createOperationLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::OpenCL));
    this->OpLTwoInner.reset(sgpp::op_factory::createOperationLTwoDotProductVectorized(
        *this->InnerGrid, sgpp::parallel::OpenCL));
    this->OpLTwoBound.reset(sgpp::op_factory::createOperationLTwoDotProductVectorized(
        *this->BoundGrid, sgpp::parallel::OpenCL));
    this->OpLTwoDotLaplaceInner.reset(sgpp::op_factory::createOperationLTwoDotLaplaceVectorized(
        *this->InnerGrid, sgpp::parallel::OpenCL));
    this->OpLTwoDotLaplaceBound.reset(sgpp::op_factory::createOperationLTwoDotLaplaceVectorized(
        *this->BoundGrid, sgpp::parallel::OpenCL));
#endif
  } else {
    throw sgpp::base::algorithm_exception(
        "PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::"
        "PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization "
        "was selected!");
  }

  // right hand side if System
  this->rhs = new sgpp::base::DataVector(1);
}

HeatEquationParabolicPDESolverSystemVectorizedMPI::
    ~HeatEquationParabolicPDESolverSystemVectorizedMPI() {
  delete this->BoundaryUpdate;
  delete this->GridConverter;

  if (this->InnerGrid != NULL) {
    delete this->InnerGrid;
  }

  if (this->alpha_inner != NULL) {
    delete this->alpha_inner;
  }

  delete this->rhs;
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);
  // Apply the Laplace operator
  this->OpLaplaceBound->mult(alpha, temp);
  result.axpy(-0.5, temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());
  result.setAll(0.0);

  // Apply the Laplace operator
  this->OpLaplaceInner->mult(alpha, temp);
  result.axpy(-0.5, temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixComplete(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoBound->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  // Apply the mass matrix
  this->OpLTwoInner->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorInner(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  this->OpLTwoDotLaplaceInner->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorBound(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(alpha.getSize());

  result.setAll(0.0);

  this->OpLTwoDotLaplaceBound->mult(alpha, temp);

  result.add(temp);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::setTimestepCoefficientInner(
    double timestep_coefficient) {
  this->OpLTwoDotLaplaceInner->setTimestepCoeff(timestep_coefficient);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::setTimestepCoefficientBound(
    double timestep_coefficient) {
  this->OpLTwoDotLaplaceBound->setTimestepCoeff(timestep_coefficient);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::finishTimestep() {
  // Replace the inner coefficients on the boundary grid
  this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::coarsenAndRefine(bool isLastTimestep) {}

void HeatEquationParabolicPDESolverSystemVectorizedMPI::startTimestep() {}
}  // namespace parallel
}  // namespace sgpp
