// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLogLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLogLinear::OperationGammaLogLinear(sgpp::base::GridStorage* storage,
                                                 sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinear::~OperationGammaLogLinear() {}

void OperationGammaLogLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                 size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                   size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOne(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimTwo(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  result.setAll(0.0);
}

void OperationGammaLogLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim) {
  sgpp::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
