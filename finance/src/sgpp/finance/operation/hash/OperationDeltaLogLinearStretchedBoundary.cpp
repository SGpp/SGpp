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


namespace SGPP {
namespace finance {

OperationDeltaLogLinearStretchedBoundary::OperationDeltaLogLinearStretchedBoundary(
  SGPP::base::GridStorage* storage,
  SGPP::base::DataVector& coef) : SGPP::pde::UpDownOneOpDim(storage, coef) {
}

OperationDeltaLogLinearStretchedBoundary::~OperationDeltaLogLinearStretchedBoundary() {
}

void OperationDeltaLogLinearStretchedBoundary::up(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretchedBoundary> s(func,
      this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::down(SGPP::base::DataVector&
    alpha, SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  SGPP::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretchedBoundary> s(func,
      this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::upOpDim(SGPP::base::DataVector&
    alpha, SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLogLinearStretchedBoundary::downOpDim(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

}
}