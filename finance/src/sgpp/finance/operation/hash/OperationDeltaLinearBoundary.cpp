// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationDeltaLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationDeltaLinearBoundary::OperationDeltaLinearBoundary(sgpp::base::GridStorage* storage,
                                                           sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationDeltaLinearBoundary::~OperationDeltaLinearBoundary() {}

void OperationDeltaLinearBoundary::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::down(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::upOpDim(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::downOpDim(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
