// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/OperationEllipticPDESolverSystemDirichlet.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sstream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationEllipticPDESolverSystemDirichlet::OperationEllipticPDESolverSystemDirichlet(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs) {
      this->BoundaryUpdate = new SGPP::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new SGPP::base::DirichletGridConverter();

      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->rhs, &(this->InnerGrid), &(this->rhs_inner));

      this->numGridpointsInner = this->InnerGrid->getSize();

      this->alpha_inner = NULL;
    }

    OperationEllipticPDESolverSystemDirichlet::~OperationEllipticPDESolverSystemDirichlet() {
      delete this->alpha_inner;
      delete this->rhs_inner;
      delete this->InnerGrid;
      delete this->BoundaryUpdate;
      delete this->GridConverter;
    }

    void OperationEllipticPDESolverSystemDirichlet::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      applyLOperatorInner(alpha, result);
    }

    SGPP::base::DataVector* OperationEllipticPDESolverSystemDirichlet::generateRHS() {
      if (this->InnerGrid != NULL) {
        SGPP::base::DataVector alpha_tmp_complete(*(this->rhs));
        SGPP::base::DataVector rhs_tmp_complete(*(this->rhs));

        this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
        applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

        this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
        this->rhs_inner->mult(-1.0);
      } else {
        throw new SGPP::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
      }

      return this->rhs_inner;
    }

    SGPP::base::DataVector* OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG() {
      if (this->InnerGrid != NULL) {
        if (this->alpha_inner != NULL) {
          delete this->alpha_inner;
        }

        this->alpha_inner = new SGPP::base::DataVector(this->InnerGrid->getSize());
        this->alpha_inner->setAll(0.0);
      } else {
        throw new SGPP::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG : No inner grid exists!");
      }

      return this->alpha_inner;
    }

    void OperationEllipticPDESolverSystemDirichlet::getSolutionBoundGrid(SGPP::base::DataVector& Solution, SGPP::base::DataVector& SolutionInner) {
      Solution = *(this->rhs);
      this->GridConverter->updateBoundaryCoefs(Solution, SolutionInner);
    }
  }
}