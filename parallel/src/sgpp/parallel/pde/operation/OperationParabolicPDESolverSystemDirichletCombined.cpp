// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/pde/operation/OperationParabolicPDESolverSystemDirichletCombined.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

OperationParabolicPDESolverSystemDirichletCombined::
    OperationParabolicPDESolverSystemDirichletCombined() {
  this->numSumGridpointsInner = 0;
  this->numSumGridpointsComplete = 0;
}

OperationParabolicPDESolverSystemDirichletCombined::
    ~OperationParabolicPDESolverSystemDirichletCombined() {}

void OperationParabolicPDESolverSystemDirichletCombined::mult(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result) {
  result.setAll(0.0);

  if (this->tOperationMode == "ImEul") {
    result.setAll(0.0);

    // Combined
    sgpp::base::DataVector temp(result.getSize());

    // this->OpLTwoDotLaplaceInner->setTimestepCoeff(-0.5*(-1.0)*this->TimestepSize);
    setTimestepCoefficientInner(-0.5 * (-1.0) * this->TimestepSize);

    applyMassMatrixLOperatorInner(alpha, temp);
    result.add(temp);

  } else if (this->tOperationMode == "CrNic") {
    result.setAll(0.0);
    sgpp::base::DataVector temp(result.getSize());

    // this->OpLTwoDotLaplaceInner->setTimestepCoeff(-0.5*(-0.5)*this->TimestepSize);
    setTimestepCoefficientInner(-0.5 * (-0.5) * this->TimestepSize);

    applyMassMatrixLOperatorInner(alpha, temp);

    result.add(temp);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::mult : An unknown operation mode was specified!");
  }
}

sgpp::base::DataVector* OperationParabolicPDESolverSystemDirichletCombined::generateRHS() {
  sgpp::base::DataVector rhs_complete(this->alpha_complete->getSize());

  if (this->tOperationMode == "ImEul") {
    rhs_complete.setAll(0.0);

    applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
  } else if (this->tOperationMode == "CrNic") {
    rhs_complete.setAll(0.0);

    sgpp::base::DataVector temp(rhs_complete.getSize());
    sgpp::base::DataVector temp2(rhs_complete.getSize());
    sgpp::base::DataVector myAlpha(*this->alpha_complete);

    // this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(0.5)*this->TimestepSize);
    setTimestepCoefficientBound(-0.5 * (0.5) * this->TimestepSize);

    applyMassMatrixLOperatorBound(myAlpha, temp);

    rhs_complete.add(temp);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was "
        "specified!");
  }

  // Now we have the right hand side, lets apply the riskfree rate for the next timestep
  this->startTimestep();

  // Now apply the boundary ansatzfunctions to the inner ansatzfunctions
  sgpp::base::DataVector result_complete(this->alpha_complete->getSize());
  sgpp::base::DataVector alpha_bound(*this->alpha_complete);

  result_complete.setAll(0.0);

  this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

  // apply CG Matrix
  if (this->tOperationMode == "ImEul") {
    sgpp::base::DataVector temp(alpha_bound.getSize());
    sgpp::base::DataVector temp2(alpha_bound.getSize());

    // this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(-1.0)*this->TimestepSize);
    setTimestepCoefficientBound(-0.5 * (-1.0) * this->TimestepSize);

    applyMassMatrixLOperatorBound(alpha_bound, temp);

    result_complete.add(temp);
  } else if (this->tOperationMode == "CrNic") {
    sgpp::base::DataVector temp(alpha_bound.getSize());
    sgpp::base::DataVector temp2(alpha_bound.getSize());

    // this->OpLTwoDotLaplaceBound->setTimestepCoeff(-0.5*(-0.5)*this->TimestepSize);
    setTimestepCoefficientBound(-0.5 * (-0.5) * this->TimestepSize);

    applyMassMatrixLOperatorBound(alpha_bound, temp);

    result_complete.add(temp);
  } else {
    throw sgpp::base::algorithm_exception(
        "OperationParabolicPDESolverSystem::generateRHS : An unknown operation mode was "
        "specified!");
  }

  rhs_complete.sub(result_complete);

  if (this->rhs != NULL) {
    delete this->rhs;
  }

  this->rhs = new sgpp::base::DataVector(this->alpha_inner->getSize());
  this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

  return this->rhs;
}

sgpp::base::DataVector*
OperationParabolicPDESolverSystemDirichletCombined::getGridCoefficientsForCG() {
  this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

  return this->alpha_inner;
}
}  // namespace parallel
}  // namespace sgpp
