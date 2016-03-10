// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonFLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationHestonFLinear::OperationHestonFLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationHestonFLinear::~OperationHestonFLinear() {}

void OperationHestonFLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonFLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonFLinear::upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonFLinear::downOpDim(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
