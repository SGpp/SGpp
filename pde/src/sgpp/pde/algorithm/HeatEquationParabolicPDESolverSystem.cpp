/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

using namespace SGPP::op_factory;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    HeatEquationParabolicPDESolverSystem::HeatEquationParabolicPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, double a, double TimestepSize, std::string OperationMode) {
      this->a = a;
      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;
      this->InnerGrid = NULL;
      this->alpha_inner = NULL;

      this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new SGPP::base::DirichletGridConverter();

      this->OpLaplaceBound = createOperationLaplace(SparseGrid);
      this->OpMassBound = SGPP::op_factory::createOperationLTwoDotProduct(SparseGrid);

      // create the inner grid
      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

      //Create needed operations, on inner grid
      this->OpLaplaceInner = createOperationLaplace(*this->InnerGrid);
      this->OpMassInner = SGPP::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);

      // right hand side if System
      this->rhs = new SGPP::base::DataVector(1);
    }

    HeatEquationParabolicPDESolverSystem::~HeatEquationParabolicPDESolverSystem() {
      delete this->OpLaplaceBound;
      delete this->OpMassBound;
      delete this->OpLaplaceInner;
      delete this->OpMassInner;

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

    void HeatEquationParabolicPDESolverSystem::applyMassMatrixComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);

      SGPP::base::DataVector temp(alpha.getSize());

      // Apply the mass matrix
      this->OpMassBound->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystem::applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);

      SGPP::base::DataVector temp(alpha.getSize());

      // Apply the laplace Operator rate
      this->OpLaplaceBound->mult(alpha, temp);
      result.axpy((-1.0)*this->a, temp);
    }

    void HeatEquationParabolicPDESolverSystem::applyMassMatrixInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);

      SGPP::base::DataVector temp(alpha.getSize());

      // Apply the mass matrix
      this->OpMassInner->mult(alpha, temp);

      result.add(temp);
    }

    void HeatEquationParabolicPDESolverSystem::applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      result.setAll(0.0);

      SGPP::base::DataVector temp(alpha.getSize());

      // Apply the laplace Operator rate
      this->OpLaplaceInner->mult(alpha, temp);
      result.axpy((-1.0)*this->a, temp);
    }

    void HeatEquationParabolicPDESolverSystem::finishTimestep() {
      // Replace the inner coefficients on the boundary grid
      this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
    }

    void HeatEquationParabolicPDESolverSystem::coarsenAndRefine(bool isLastTimestep) {
    }

    void HeatEquationParabolicPDESolverSystem::startTimestep() {
    }

  }
}
