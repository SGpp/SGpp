// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationDeltaLogLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationDeltaLogLinearStretched::OperationDeltaLogLinearStretched(sgpp::base::GridStorage* storage,
                                                                   sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationDeltaLogLinearStretched::~OperationDeltaLogLinearStretched() {}

void OperationDeltaLogLinearStretched::up(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::down(sgpp::base::DataVector& alpha,
                                            sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::upOpDim(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::downOpDim(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
