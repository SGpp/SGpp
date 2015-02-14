// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationParabolicPDESolverSystemFreeBoundaries::OperationParabolicPDESolverSystemFreeBoundaries() {
      this->numSumGridpointsInner = 0;
      this->numSumGridpointsComplete = 0;
    }

    OperationParabolicPDESolverSystemFreeBoundaries::~OperationParabolicPDESolverSystemFreeBoundaries() {
    }

    void OperationParabolicPDESolverSystemFreeBoundaries::mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      if (this->tOperationMode == "ExEul") {
        applyMassMatrix(alpha, result);
      } else if (this->tOperationMode == "ImEul") {
        result.setAll(0.0);

        SGPP::base::DataVector temp(alpha.getSize());

        applyMassMatrix(alpha, temp);
        result.add(temp);

        applyLOperator(alpha, temp);
        result.axpy((-1.0)*this->TimestepSize, temp);
      } else if (this->tOperationMode == "CrNic") {
        result.setAll(0.0);

        SGPP::base::DataVector temp(alpha.getSize());

        applyMassMatrix(alpha, temp);
        result.add(temp);

        applyLOperator(alpha, temp);
        result.axpy((-0.5)*this->TimestepSize, temp);
      } else if (this->tOperationMode == "AdBas") {
        result.setAll(0.0);

        applyMassMatrix(alpha, result);
      } else if (this->tOperationMode == "BDF2") {
        float_t tDiff = this->TimestepSize / this->TimestepSize_old;
        float_t alpha0 = (2.0 * tDiff + 1.0) / (tDiff + 1.0);
        result.setAll(0.0);

        SGPP::base::DataVector temp(alpha.getSize());

        applyMassMatrix(alpha, temp);

        temp.mult(alpha0);
        result.add(temp);

        applyLOperator(alpha, temp);
        result.axpy((-1.0)*this->TimestepSize, temp);
      } else if (this->tOperationMode == "F23") {
        result.setAll(0.0);
        float_t tDiff = this->TimestepSize / this->TimestepSize_old;
        float_t alpha0 = 1.0 / (1.0 + tDiff);

        applyMassMatrix(alpha, result);
        result.mult(alpha0);

      } else {
        throw new SGPP::base::algorithm_exception("OperationParabolicPDESolverSystemNeumann::mult : An unknown operation mode was specified!");
      }
    }

    SGPP::base::DataVector* OperationParabolicPDESolverSystemFreeBoundaries::generateRHS() {
      SGPP::base::DataVector rhs_complete(this->alpha_complete->getSize());

      if (this->tOperationMode == "ExEul") {
        rhs_complete.setAll(0.0);

        SGPP::base::DataVector temp(this->alpha_complete->getSize());

        applyMassMatrix(*this->alpha_complete, temp);
        rhs_complete.add(temp);

        applyLOperator(*this->alpha_complete, temp);
        rhs_complete.axpy(this->TimestepSize, temp);
      } else if (this->tOperationMode == "ImEul") {
        rhs_complete.setAll(0.0);

        applyMassMatrix(*this->alpha_complete, rhs_complete);
      } else if (this->tOperationMode == "CrNic") {
        rhs_complete.setAll(0.0);

        SGPP::base::DataVector temp(this->alpha_complete->getSize());

        applyMassMatrix(*this->alpha_complete, temp);
        rhs_complete.add(temp);

        applyLOperator(*this->alpha_complete, temp);
        rhs_complete.axpy((0.5)*this->TimestepSize, temp);
      } else if (this->tOperationMode == "AdBas") {
        rhs_complete.setAll(0.0);

        SGPP::base::DataVector temp(this->alpha_complete->getSize());

        applyMassMatrix(*this->alpha_complete, temp);
        rhs_complete.add(temp);

        applyLOperator(*this->alpha_complete, temp);

        temp.mult((2.0) + this->TimestepSize / this->TimestepSize_old);

        SGPP::base::DataVector temp_old(this->alpha_complete->getSize());
        applyMassMatrix(*this->alpha_complete_old, temp_old);
        applyLOperator(*this->alpha_complete_old, temp_old);
        temp_old.mult(this->TimestepSize / this->TimestepSize_old);
        temp.sub(temp_old);

        rhs_complete.axpy((0.5)*this->TimestepSize, temp);
      } else if (this->tOperationMode == "BDF2") {
        rhs_complete.setAll(0.0);

        SGPP::base::DataVector temp(this->alpha_complete->getSize());

        applyMassMatrix(*this->alpha_complete, temp);

        float_t tDiff = this->TimestepSize / this->TimestepSize_old;

        float_t alpha1 = tDiff + 1.0;
        temp.mult(alpha1);
        rhs_complete.add(temp);

        SGPP::base::DataVector temp_old(this->alpha_complete->getSize());
        applyMassMatrix(*this->alpha_complete_old, temp_old);

        float_t alpha2 = tDiff * tDiff / (1.0 + tDiff);
        temp_old.mult(alpha2);
        rhs_complete.sub(temp_old);
      } else if (this->tOperationMode == "F23") {
        rhs_complete.setAll(0.0);
        float_t tDiff = this->TimestepSize / this->TimestepSize_old;
        float_t alpha0 = (1.0 + tDiff);
        float_t alpha1 = alpha0 * (tDiff - 1.0);
        float_t alpha2 = -alpha0 * (tDiff * tDiff / (tDiff + 1.0));


        SGPP::base::DataVector temp(this->alpha_complete->getSize());
        SGPP::base::DataVector temp_old(this->alpha_complete->getSize());

        applyMassMatrix(*this->alpha_complete, temp);
        temp.mult(alpha1);
        rhs_complete.sub(temp);

        applyLOperator(*this->alpha_complete, temp);
        temp.mult(alpha0 * this->TimestepSize);

        applyMassMatrix(*this->alpha_complete_old, temp_old);
        temp_old.mult(alpha2);
        rhs_complete.sub(temp_old);

        rhs_complete.add(temp);
      } else {
        throw new SGPP::base::algorithm_exception("OperationParabolicPDESolverSystemNeumann::generateRHS : An unknown operation mode was specified!");
      }

      // Now we have the right hand side, lets apply the riskfree rate for the next timestep
      this->startTimestep();

      if (this->rhs != NULL) {
        delete this->rhs;
      }

      this->rhs = new SGPP::base::DataVector(rhs_complete);

      return this->rhs;
    }

    SGPP::base::DataVector* OperationParabolicPDESolverSystemFreeBoundaries::getGridCoefficientsForCG() {
      return this->alpha_complete;
    }

  }
}
