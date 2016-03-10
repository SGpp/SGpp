// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonZLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationHestonZLinear::OperationHestonZLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationHestonZLinear::~OperationHestonZLinear() {}

void OperationHestonZLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::downOpDim(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
