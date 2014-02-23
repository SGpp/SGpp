/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemVectorizedMPI.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"

#include "base/exception/algorithm_exception.hpp"

#include <cstring>

namespace sg {
  namespace parallel {

    HeatEquationParabolicPDESolverSystemVectorizedMPI::HeatEquationParabolicPDESolverSystemVectorizedMPI(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode) {
      this->a = a;
      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;
      this->InnerGrid = NULL;
      this->alpha_inner = NULL;

      this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new sg::base::DirichletGridConverter();

      // create the inner grid
      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

      // Create operations
      char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

      if (! strcmp(alg_selector, "X86SIMD")) {
        this->OpLaplaceInner = sg::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, sg::parallel::X86SIMD);
        this->OpLaplaceBound = sg::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, sg::parallel::X86SIMD);
        this->OpLTwoInner = sg::op_factory::createOperationLTwoDotProductVectorized(*this->InnerGrid, sg::parallel::X86SIMD);
        this->OpLTwoBound = sg::op_factory::createOperationLTwoDotProductVectorized(*this->BoundGrid, sg::parallel::X86SIMD);
        this->OpLTwoDotLaplaceInner = sg::op_factory::createOperationLTwoDotLaplaceVectorized(*this->InnerGrid, sg::parallel::X86SIMD);
        this->OpLTwoDotLaplaceBound = sg::op_factory::createOperationLTwoDotLaplaceVectorized(*this->BoundGrid, sg::parallel::X86SIMD);
#ifdef USEOCL
      } else if (! strcmp(alg_selector, "OCL")) {
        this->OpLaplaceInner = sg::op_factory::createOperationLaplaceVectorized(*this->InnerGrid, sg::parallel::OpenCL);
        this->OpLaplaceBound = sg::op_factory::createOperationLaplaceVectorized(*this->BoundGrid, sg::parallel::OpenCL);
        this->OpLTwoInner = sg::op_factory::createOperationLTwoDotProductVectorized(*this->InnerGrid, sg::parallel::OpenCL);
        this->OpLTwoBound = sg::op_factory::createOperationLTwoDotProductVectorized(*this->BoundGrid, sg::parallel::OpenCL);
        this->OpLTwoDotLaplaceInner = sg::op_factory::createOperationLTwoDotLaplaceVectorized(*this->InnerGrid, sg::parallel::OpenCL);
        this->OpLTwoDotLaplaceBound = sg::op_factory::createOperationLTwoDotLaplaceVectorized(*this->BoundGrid, sg::parallel::OpenCL);
#endif
      } else {
        throw sg::base::algorithm_exception("PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI::PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : no supported vectorization was selected!");
      }

      // right hand side if System
      this->rhs = new sg::base::DataVector(1);
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

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);
      // Apply the Laplace operator
      this->OpLaplaceBound->mult(alpha, temp);
      result.axpy(-0.5, temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());
      result.setAll(0.0);

      // Apply the Laplace operator
      this->OpLaplaceInner->mult(alpha, temp);
      result.axpy(-0.5, temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoBound->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoInner->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      this->OpLTwoDotLaplaceInner->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystemVectorizedMPI::applyMassMatrixLOperatorBound(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

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
