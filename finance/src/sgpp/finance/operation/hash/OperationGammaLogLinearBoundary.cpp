// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLogLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace SGPP {
namespace finance {

OperationGammaLogLinearBoundary::OperationGammaLogLinearBoundary(SGPP::base::GridStorage* storage,
                                                                 SGPP::base::DataMatrix& coef)
    : SGPP::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinearBoundary::~OperationGammaLogLinearBoundary() {}

void OperationGammaLogLinearBoundary::up(SGPP::base::DataVector& alpha,
                                         SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::down(SGPP::base::DataVector& alpha,
                                           SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimOne(SGPP::base::DataVector& alpha,
                                                 SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<PhidPhiUpBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimOne(SGPP::base::DataVector& alpha,
                                                   SGPP::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<PhidPhiDownBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimTwo(SGPP::base::DataVector& alpha,
                                                 SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimTwo(SGPP::base::DataVector& alpha,
                                                   SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                            SGPP::base::DataVector& result,
                                                            size_t dim) {
  SGPP::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                                              SGPP::base::DataVector& result,
                                                              size_t dim) {
  SGPP::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
