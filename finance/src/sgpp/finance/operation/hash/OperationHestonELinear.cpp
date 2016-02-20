// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonELinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

OperationHestonELinear::OperationHestonELinear(SGPP::base::GridStorage* storage,
                                               SGPP::base::DataVector& coef)
    : SGPP::pde::UpDownOneOpDim(storage, coef) {}

OperationHestonELinear::~OperationHestonELinear() {}

void OperationHestonELinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonELinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonELinear::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                     size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonELinear::downOpDim(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
