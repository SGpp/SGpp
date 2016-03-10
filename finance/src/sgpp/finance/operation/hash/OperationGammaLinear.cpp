// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLinear::OperationGammaLinear(sgpp::base::GridStorage* storage,
                                           sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLinear::~OperationGammaLinear() {}

void OperationGammaLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                              size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOne(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XPhidPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XPhidPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
