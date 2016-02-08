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

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

OperationGammaLinear::OperationGammaLinear(SGPP::base::GridStorage* storage,
    SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
}

OperationGammaLinear::~OperationGammaLinear() {
}

void OperationGammaLinear::up(SGPP::base::DataVector& alpha,
                              SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::down(SGPP::base::DataVector& alpha,
                                SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOne(SGPP::base::DataVector& alpha,
                                      SGPP::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<XPhidPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOne(SGPP::base::DataVector& alpha,
                                        SGPP::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<XPhidPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimTwo(SGPP::base::DataVector& alpha,
                                      SGPP::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimTwo(SGPP::base::DataVector& alpha,
                                        SGPP::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<SqXdPhidPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLinear::downOpDimOneAndOpDimTwo(SGPP::base::DataVector&
    alpha, SGPP::base::DataVector& result, size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<SqXdPhidPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

}
}