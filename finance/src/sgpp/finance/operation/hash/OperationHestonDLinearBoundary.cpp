// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonDLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationHestonDLinearBoundary::OperationHestonDLinearBoundary(sgpp::base::GridStorage* storage,
                                                               sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationHestonDLinearBoundary::~OperationHestonDLinearBoundary() {}

void OperationHestonDLinearBoundary::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonDLinearBoundary::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonDLinearBoundary::upOpDim(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * dphi
  XdPhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationHestonDLinearBoundary::downOpDim(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * dphi
  XdPhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
