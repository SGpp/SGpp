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

namespace SGPP {
namespace finance {

OperationGammaLogLinearStretched::OperationGammaLogLinearStretched(SGPP::base::GridStorage* storage,
                                                                   SGPP::base::DataMatrix& coef)
    : SGPP::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinearStretched::~OperationGammaLogLinearStretched() {}

void OperationGammaLogLinearStretched::up(SGPP::base::DataVector& alpha,
                                          SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::down(SGPP::base::DataVector& alpha,
                                            SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOne(SGPP::base::DataVector& alpha,
                                                  SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<PhidPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimOne(SGPP::base::DataVector& alpha,
                                                    SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<PhidPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimTwo(SGPP::base::DataVector& alpha,
                                                  SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretched func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::downOpDimTwo(SGPP::base::DataVector& alpha,
                                                    SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretched func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationGammaLogLinearStretched::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                             SGPP::base::DataVector& result,
                                                             size_t dim) {
  result.setAll(0.0);
}

void OperationGammaLogLinearStretched::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                               SGPP::base::DataVector& result,
                                                               size_t dim) {
  SGPP::pde::DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
