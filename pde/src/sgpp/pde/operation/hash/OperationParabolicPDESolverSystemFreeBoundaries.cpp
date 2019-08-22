// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationParabolicPDESolverSystemFreeBoundaries::OperationParabolicPDESolverSystemFreeBoundaries() {
  this->numSumGridpointsInner = 0;
  this->numSumGridpointsComplete = 0;
}

OperationParabolicPDESolverSystemFreeBoundaries::
    ~OperationParabolicPDESolverSystemFreeBoundaries() {}

void OperationParabolicPDESolverSystemFreeBoundaries::mult(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result) {
  if (this->tOperationMode == "ExEul") {
    applyMassMatrix(alpha, result);
  } else if (this->tOperationMode == "ImEul") {
    result.setAll(0.0);

    sgpp::base::DataVector temp(alpha.getSize());

    applyMassMatrix(alpha, temp);
    result.add(temp);

    applyLOperator(alpha, temp);
    result.axpy((-1.0) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "CrNic") {
    result.setAll(0.0);

    sgpp::base::DataVector temp(alpha.getSize());

    applyMassMatrix(alpha, temp);
    result.add(temp);

    applyLOperator(alpha, temp);
    result.axpy((-0.5) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "AdBas") {
    result.setAll(0.0);

    applyMassMatrix(alpha, result);
  } else if (this->tOperationMode == "BDF2") {
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = (2.0 * tDiff + 1.0) / (tDiff + 1.0);
    result.setAll(0.0);

    sgpp::base::DataVector temp(alpha.getSize());

    applyMassMatrix(alpha, temp);

    temp.mult(alpha0);
    result.add(temp);

    applyLOperator(alpha, temp);
    result.axpy((-1.0) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "F23") {
    result.setAll(0.0);
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = 1.0 / (1.0 + tDiff);

    applyMassMatrix(alpha, result);
    result.mult(alpha0);

  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystemNeumann::mult : An unknown operation mode was "
        "specified!");
  }
}

sgpp::base::DataVector* OperationParabolicPDESolverSystemFreeBoundaries::generateRHS() {
  sgpp::base::DataVector rhs_complete(this->alpha_complete->getSize());

  if (this->tOperationMode == "ExEul") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrix(*this->alpha_complete, temp);
    rhs_complete.add(temp);

    applyLOperator(*this->alpha_complete, temp);
    rhs_complete.axpy(this->TimestepSize, temp);
  } else if (this->tOperationMode == "ImEul") {
    rhs_complete.setAll(0.0);

    applyMassMatrix(*this->alpha_complete, rhs_complete);
  } else if (this->tOperationMode == "CrNic") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrix(*this->alpha_complete, temp);
    rhs_complete.add(temp);

    applyLOperator(*this->alpha_complete, temp);
    rhs_complete.axpy((0.5) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "AdBas") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrix(*this->alpha_complete, temp);
    rhs_complete.add(temp);

    applyLOperator(*this->alpha_complete, temp);

    temp.mult((2.0) + this->TimestepSize / this->TimestepSize_old);

    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());
    applyMassMatrix(*this->alpha_complete_old, temp_old);
    applyLOperator(*this->alpha_complete_old, temp_old);
    temp_old.mult(this->TimestepSize / this->TimestepSize_old);
    temp.sub(temp_old);

    rhs_complete.axpy((0.5) * this->TimestepSize, temp);
  } else if (this->tOperationMode == "BDF2") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(this->alpha_complete->getSize());

    applyMassMatrix(*this->alpha_complete, temp);

    double tDiff = this->TimestepSize / this->TimestepSize_old;

    double alpha1 = tDiff + 1.0;
    temp.mult(alpha1);
    rhs_complete.add(temp);

    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());
    applyMassMatrix(*this->alpha_complete_old, temp_old);

    double alpha2 = tDiff * tDiff / (1.0 + tDiff);
    temp_old.mult(alpha2);
    rhs_complete.sub(temp_old);
  } else if (this->tOperationMode == "F23") {
    rhs_complete.setAll(0.0);
    double tDiff = this->TimestepSize / this->TimestepSize_old;
    double alpha0 = (1.0 + tDiff);
    double alpha1 = alpha0 * (tDiff - 1.0);
    double alpha2 = -alpha0 * (tDiff * tDiff / (tDiff + 1.0));

    sgpp::base::DataVector temp(this->alpha_complete->getSize());
    sgpp::base::DataVector temp_old(this->alpha_complete->getSize());

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
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystemNeumann::generateRHS : An unknown operation mode was "
        "specified!");
  }

  // Now we have the right hand side, lets apply the riskfree rate for the next timestep
  this->startTimestep();

  if (this->rhs != nullptr) {
    delete this->rhs;
  }

  this->rhs = new sgpp::base::DataVector(rhs_complete);

  return this->rhs;
}

sgpp::base::DataVector*
OperationParabolicPDESolverSystemFreeBoundaries::getGridCoefficientsForCG() {
  return this->alpha_complete;
}
}  // namespace pde
}  // namespace sgpp
