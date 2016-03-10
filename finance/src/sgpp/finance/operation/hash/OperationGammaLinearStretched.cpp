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

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLinearStretched::OperationGammaLinearStretched(sgpp::base::GridStorage* storage,
                                                             sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLinearStretched::~OperationGammaLinearStretched() {}

void OperationGammaLinearStretched::up(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::down(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOne(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<XPhidPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOne(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<XPhidPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimTwo(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                          sgpp::base::DataVector& result,
                                                          size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinearStretched::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
