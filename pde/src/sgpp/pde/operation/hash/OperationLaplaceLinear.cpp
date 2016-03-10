// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceLinear::OperationLaplaceLinear(sgpp::base::GridStorage* storage)
    : UpDownOneOpDim(storage) {}

OperationLaplaceLinear::OperationLaplaceLinear(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataVector& coef)
    : UpDownOneOpDim(storage, coef) {}

OperationLaplaceLinear::~OperationLaplaceLinear() {}

void OperationLaplaceLinear::specialOP(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim,
                                       size_t gradient_dim) {
  // In direction gradient_dim we only calculate the norm of the gradient
  // The up-part is empty, thus omitted
  if (dim > 0) {
    sgpp::base::DataVector temp(alpha.getSize());
    updown(alpha, temp, dim - 1, gradient_dim);
    downOpDim(temp, result, gradient_dim);
  } else {
    // Terminates dimension recursion
    downOpDim(alpha, result, gradient_dim);
  }
}

void OperationLaplaceLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                size_t dim) {
  PhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                  size_t dim) {
  PhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::downOpDim(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) {
  DowndPhidPhiBBIterativeLinear myDown(this->storage);
  myDown(alpha, result, dim);
}

void OperationLaplaceLinear::upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {}
}  // namespace pde
}  // namespace sgpp
