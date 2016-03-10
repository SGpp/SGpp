// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLogLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/PhidPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/PhidPhiUpBBLinearStretched.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp>
#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiUpBBLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/DowndPhidPhiBBIterativeLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLogLinearStretched::OperationGammaLogLinearStretched(sgpp::base::GridStorage* storage,
                                                                   sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinearStretched::~OperationGammaLogLinearStretched() {}

void OperationGammaLogLinearStretched::up(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::down(sgpp::base::DataVector& alpha,
                                            sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOne(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimOne(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result,
                                                             size_t dim) {
  result.setAll(0.0);
}

void OperationGammaLogLinearStretched::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result,
                                                               size_t dim) {
  sgpp::pde::DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
