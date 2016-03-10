// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationDeltaLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/XdPhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationDeltaLinearStretchedBoundary::OperationDeltaLinearStretchedBoundary(
    sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationDeltaLinearStretchedBoundary::~OperationDeltaLinearStretchedBoundary() {}

void OperationDeltaLinearStretchedBoundary::up(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::down(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::upOpDim(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearStretchedBoundary::downOpDim(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim) {
  // x * dphi * phi
  XdPhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
