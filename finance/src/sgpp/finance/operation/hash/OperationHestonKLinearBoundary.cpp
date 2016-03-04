// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonKLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationHestonKLinearBoundary::OperationHestonKLinearBoundary(sgpp::base::GridStorage* storage,
                                                               double***** coef)
    : sgpp::pde::UpDownFourOpDims(storage, coef) {}

OperationHestonKLinearBoundary::~OperationHestonKLinearBoundary() {}

// Unidirectional
void OperationHestonKLinearBoundary::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
void OperationHestonKLinearBoundary::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

// Singles
void OperationHestonKLinearBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::downOpDimThree(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::upOpDimThree(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::downOpDimFour(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonKLinearBoundary::upOpDimFour(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // sqrtX phi phi
  SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

// Doubles
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result,
                                                           size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result,
                                                               size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result,
                                                              size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result,
                                                               size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result,
                                                              size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                sgpp::base::DataVector& result,
                                                                size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result,
                                                              size_t dim) {}

// Triples
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThree(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThree(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}

// Quadruples
void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {}
}  // namespace finance
}  // namespace sgpp
