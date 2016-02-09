// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

OperationGammaLinearStretched::OperationGammaLinearStretched(
  SGPP::base::GridStorage* storage,
  SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
}

OperationGammaLinearStretched::~OperationGammaLinearStretched() {
}

void OperationGammaLinearStretched::up(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::down(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretched> s(func,
      this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOne(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<XPhidPhiUpBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOne(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<XPhidPhiDownBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimTwo(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimTwo(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOneAndOpDimTwo(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<SqXdPhidPhiUpBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOneAndOpDimTwo(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<SqXdPhidPhiDownBBLinearStretched> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

}
}