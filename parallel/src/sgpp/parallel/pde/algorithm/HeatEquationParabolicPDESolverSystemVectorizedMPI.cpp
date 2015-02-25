// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemVectorizedMPI.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    HeatEquationParabolicPDESolverSystemVectorizedMPI::HeatEquationParabolicPDESolverSystemVectorizedMPI(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode) {
      this->a = a;
      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;
      this->InnerGrid = NULL;
      this->alpha_inner = NULL;

      this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new SGPP::base::DirichletGridConverter();

      // create the inner grid
      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

      // Create operations
      char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

      if (! strcmp(alg_selector, "X86SIMD")) {
        this->OpLaplaceInner = SGPP::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, SGPP::parallel::X86SIMD);
        this->OpLaplaceBound = SGPP::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, SGPP::parallel::X86SIMD);
        this->OpLTwoInner = SGPP::op_factory::createOperationLTwoDotProductVectorized(*this->InnerGrid, SGPP::parallel::X86SIMD);
        this->OpLTwoBound = SGPP::op_factory::createOperationLTwoDotProductVectorized(*this->BoundGrid, SGPP::parallel::X86SIMD);
        this->OpLTwoDotLaplaceInner = SGPP::op_factory::createOperationLTwoDotLaplaceVectorized(*this->InnerGrid, SGPP::parallel::X86SIMD);
        this->OpLTwoDotLaplaceBound = SGPP::op_factory::createOperationLTwoDotLaplaceVectorized(*this->BoundGrid, SGPP::parallel::X86SIMD);
#ifdef USEOCL
      } else if (! strcmp(alg_selector, "OCL")) {
        this->OpLaplaceInner = SGPP::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, SGPP::parallel::OpenCL);
        this->OpLaplaceBound = SGPP::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, SGPP::parallel::OpenCL);
        this->OpLTwoInner = SGPP::op_factory::createOperationLTwoDotProductVectorized(*this->InnerGrid, SGPP::parallel::OpenCL);
        this->OpLTwoBound = SGPP::op_factory::createOperationLTwoDotProductVectorized(*this->BoundGrid, SGPP::parallel::OpenCL);
        this->OpLTwoDotLaplaceInner = SGPP::op_factory::createOperationLTwoDotLaplaceVectorized(*this->InnerGrid, SGPP::parallel::OpenCL);
        this->OpLTwoDotLaplaceBound = SGPP::op_factory::createOperationLTwoDotLaplaceVectorized(*this->BoundGrid, SGPP::parallel::OpenCL);
#endif
      } else {
        throw SGPP::base::algorithm_exception("PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization was selected!");
      }

      // right hand side if System
      this->rhs = new SGPP::base::DataVector(1);
    }

    HeatEquationParabolicPDESolverSystemVectorizedMPI::~HeatEquationParabolicPDESolverSystemVectorizedMPI() {
      delete this->OpLaplaceInner;
      delete this->OpLaplaceBound;
      delete this->OpLTwoInner;
      delete this->OpLTwoBound;
      delete this->OpLTwoDotLaplaceInner;
      delete this->OpLTwoDotLaplaceBound;

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

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);
      // Apply the Laplace operator
      this->OpLaplaceBound->mult(alpha, temp);
      result.axpy(-0.5, temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());
      result.setAll(0.0);

      // Apply the Laplace operator
      this->OpLaplaceInner->mult(alpha, temp);
      result.axpy(-0.5, temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoBound->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoInner->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      this->OpLTwoDotLaplaceInner->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorBound(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      SGPP::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      this->OpLTwoDotLaplaceBound->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::setTimestepCoefficientInner(double timestep_coefficient) {
      this->OpLTwoDotLaplaceInner->setTimestepCoeff(timestep_coefficient);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::setTimestepCoefficientBound(double timestep_coefficient) {
      this->OpLTwoDotLaplaceBound->setTimestepCoeff(timestep_coefficient);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::finishTimestep() {
      // Replace the inner coefficients on the boundary grid
      this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::coarsenAndRefine(bool isLastTimestep) {
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::startTimestep() {
    }

  }
}