// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/DowndPhidPhiBBIterativeLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceLinearStretched::OperationLaplaceLinearStretched(sgpp::base::GridStorage* storage)
    : UpDownOneOpDim(storage) {}

OperationLaplaceLinearStretched::~OperationLaplaceLinearStretched() {}

void OperationLaplaceLinearStretched::specialOP(sgpp::base::DataVector& alpha,
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

void OperationLaplaceLinearStretched::up(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  PhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinearStretched> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinearStretched::down(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  PhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinearStretched> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinearStretched::downOpDim(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  DowndPhidPhiBBIterativeLinearStretched myDown(this->storage);
  myDown(alpha, result, dim);
}

void OperationLaplaceLinearStretched::upOpDim(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result, size_t dim) {}
}  // namespace pde
}  // namespace sgpp
