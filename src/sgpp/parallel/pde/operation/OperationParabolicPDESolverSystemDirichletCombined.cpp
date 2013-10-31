/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#include "parallel/pde/operation/OperationParabolicPDESolverSystemDirichletCombined.hpp"
#include "base/exception/algorithm_exception.hpp"

namespace sg {
  namespace parallel {

    OperationParabolicPDESolverSystemDirichletCombined::OperationParabolicPDESolverSystemDirichletCombined() {
      this->numSumGridpointsInner = 0;
      this->numSumGridpointsComplete = 0;
    }

    OperationParabolicPDESolverSystemDirichletCombined::~OperationParabolicPDESolverSystemDirichletCombined() {
    }

    void OperationParabolicPDESolverSystemDirichletCombined::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      if (this->tOperationMode == "ImEul") {
        result.setAll(0.0);

        // Combined
        sg::base::DataVector temp(result.getSize());

        //this->OpLTwoDotLaplaceInner->setTimestepCoeff(-0.5*(-1.0)*this->TimestepSize);
        setTimestepCoefficientInner(-0.5 * (-1.0)*this->TimestepSize);

        applyMassMatrixLOperatorInner(alpha, temp);
        result.add(temp);

      } else if (this->tOperationMode == "CrNic") {
        result.setAll(0.0);
        sg::base::DataVector temp(result.getSize());

        //this->OpLTwoDotLaplaceInner->setTimestepCoeff(-0.5*(-0.5)*this->TimestepSize);
        setTimestepCoefficientInner(-0.5 * (-0.5)*this->TimestepSize);

        applyMassMatrixLOperatorInner(alpha, temp);

        result.add(temp);
      } else {
        throw new sg::base::algorithm_exception("OperationParabolicPDESolverSystem::mult : An unknown operation mode was specified!");
      }
    }

    sg::base::DataVector* OperationParabolicPDESolverSystemDirichletCombined::generateRHS() {
      sg::base::DataVector rhs_complete(this->alpha_complete->getSize());

      if (this->tOperationMode == "ImEul") {
        rhs_complete.setAll(0.0);

        applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
      } else if (this->tOperationMode == "CrNic") {
        rhs_complete.setAll(0.0);

        sg::base::DataVector temp(rhs_complete.getSize());
        sg::base::DataVector temp2(rhs_complete.getSize());
        sg::base::DataVector myAlpha(*this->alpha_complete);

        //this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(0.5)*this->TimestepSize);
        setTimestepCoefficientBound(-0.5 * (0.5)*this->TimestepSize);

        applyMassMatrixLOperatorBound(myAlpha, temp);

        rhs_complete.add(temp);
      } else {
        throw new sg::base::algorithm_exception("OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was specified!");
      }

      // Now we have the right hand side, lets apply the riskfree rate for the next timestep
      this->startTimestep();

      // Now apply the boundary ansatzfunctions to the inner ansatzfunctions
      sg::base::DataVector result_complete(this->alpha_complete->getSize());
      sg::base::DataVector alpha_bound(*this->alpha_complete);

      result_complete.setAll(0.0);

      this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

      // apply CG Matrix
      if (this->tOperationMode == "ImEul") {
        sg::base::DataVector temp(alpha_bound.getSize());
        sg::base::DataVector temp2(alpha_bound.getSize());

        //this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(-1.0)*this->TimestepSize);
        setTimestepCoefficientBound(-0.5 * (-1.0)*this->TimestepSize);

        applyMassMatrixLOperatorBound(alpha_bound, temp);

        result_complete.add(temp);
      } else if (this->tOperationMode == "CrNic") {
        sg::base::DataVector temp(alpha_bound.getSize());
        sg::base::DataVector temp2(alpha_bound.getSize());

        //this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(-0.5)*this->TimestepSize);
        setTimestepCoefficientBound(-0.5 * (-0.5)*this->TimestepSize);

        applyMassMatrixLOperatorBound(alpha_bound, temp);

        result_complete.add(temp);
      } else {
        throw new sg::base::algorithm_exception("OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was specified!");
      }

      rhs_complete.sub(result_complete);

      if (this->rhs != NULL) {
        delete this->rhs;
      }

      this->rhs = new sg::base::DataVector(this->alpha_inner->getSize());
      this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

      return this->rhs;
    }

    sg::base::DataVector* OperationParabolicPDESolverSystemDirichletCombined::getGridCoefficientsForCG() {
      this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

      return this->alpha_inner;
    }
  }
}
