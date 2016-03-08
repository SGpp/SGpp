// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XPhidPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XPhidPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLinearStretchedBoundary::OperationGammaLinearStretchedBoundary(
    sgpp::base::GridStorage* storage, sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLinearStretchedBoundary::~OperationGammaLinearStretchedBoundary() {}

void OperationGammaLinearStretchedBoundary::up(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::down(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
  // x * phi * dphi
  XPhidPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                         sgpp::base::DataVector& result,
                                                         size_t dim) {
  // x * phi * dphi
  XPhidPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XPhidPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                         sgpp::base::DataVector& result,
                                                         size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                                  sgpp::base::DataVector& result,
                                                                  size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearStretchedBoundary::downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                                    sgpp::base::DataVector& result,
                                                                    size_t dim) {
  // x^2 * dphi * dphi
  SqXdPhidPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<SqXdPhidPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
