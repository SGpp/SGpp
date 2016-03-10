// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonKLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonKLinear::OperationHestonKLinear(sgpp::base::GridStorage* storage, double***** coef)
    : sgpp::pde::UpDownFourOpDims(storage, coef) {}

OperationHestonKLinear::~OperationHestonKLinear() {}

// Unidirectional
void OperationHestonKLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
void OperationHestonKLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

// Singles
void OperationHestonKLinear::downOpDimOne(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::upOpDimOne(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::downOpDimTwo(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::upOpDimTwo(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::downOpDimThree(sgpp::base::DataVector& alpha,
                                            sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::upOpDimThree(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::downOpDimFour(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationHestonKLinear::upOpDimFour(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

// Doubles
void OperationHestonKLinear::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
}
void OperationHestonKLinear::upOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
}
void OperationHestonKLinear::upOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                        sgpp::base::DataVector& result,
                                                        size_t dim) {}
void OperationHestonKLinear::upOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim) {}

// Triples
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                                  sgpp::base::DataVector& result,
                                                                  size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                                sgpp::base::DataVector& result,
                                                                size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                 sgpp::base::DataVector& result,
                                                                 size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result,
                                                               size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                   sgpp::base::DataVector& result,
                                                                   size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                 sgpp::base::DataVector& result,
                                                                 size_t dim) {}
void OperationHestonKLinear::downOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                   sgpp::base::DataVector& result,
                                                                   size_t dim) {}
void OperationHestonKLinear::upOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                 sgpp::base::DataVector& result,
                                                                 size_t dim) {}

// Quadruples
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
}  // namespace finance
}  // namespace sgpp
