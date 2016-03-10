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

namespace sgpp {
namespace finance {

OperationGammaLogLinearBoundary::OperationGammaLogLinearBoundary(sgpp::base::GridStorage* storage,
                                                                 sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinearBoundary::~OperationGammaLogLinearBoundary() {}

void OperationGammaLogLinearBoundary::up(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::down(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {
  sgpp::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}

void OperationGammaLogLinearBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result,
                                                              size_t dim) {
  sgpp::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
