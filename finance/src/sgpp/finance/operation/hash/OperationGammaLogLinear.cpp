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

namespace SGPP {
namespace finance {

OperationGammaLogLinear::OperationGammaLogLinear(SGPP::base::GridStorage* storage,
                                                 SGPP::base::DataMatrix& coef)
    : SGPP::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinear::~OperationGammaLogLinear() {}

void OperationGammaLogLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                 size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                   size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOne(SGPP::base::DataVector& alpha,
                                         SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<PhidPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimOne(SGPP::base::DataVector& alpha,
                                           SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<PhidPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimTwo(SGPP::base::DataVector& alpha,
                                         SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::downOpDimTwo(SGPP::base::DataVector& alpha,
                                           SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinear::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                    SGPP::base::DataVector& result, size_t dim) {
  result.setAll(0.0);
}

void OperationGammaLogLinear::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                      SGPP::base::DataVector& result, size_t dim) {
  SGPP::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
