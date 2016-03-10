// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLinearBoundary::OperationGammaLinearBoundary(sgpp::base::GridStorage* storage,
                                                           sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLinearBoundary::~OperationGammaLinearBoundary() {}

void OperationGammaLinearBoundary::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::down(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                         sgpp::base::DataVector& result,
                                                         size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                           sgpp::base::DataVector& result,
                                                           size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
