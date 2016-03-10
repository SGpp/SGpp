// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationDeltaLogLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationDeltaLogLinearStretchedBoundary::OperationDeltaLogLinearStretchedBoundary(
    sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef)
    : sgpp::pde::UpDownOneOpDim(storage, coef) {}

OperationDeltaLogLinearStretchedBoundary::~OperationDeltaLogLinearStretchedBoundary() {}

void OperationDeltaLogLinearStretchedBoundary::up(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::down(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::upOpDim(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::downOpDim(sgpp::base::DataVector& alpha,
                                                         sgpp::base::DataVector& result,
                                                         size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
