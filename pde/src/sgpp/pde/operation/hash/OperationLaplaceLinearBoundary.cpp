// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(sgpp::base::GridStorage* storage)
    : UpDownOneOpDim(storage) {}

OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(sgpp::base::GridStorage* storage,
                                                               sgpp::base::DataVector& coef)
    : UpDownOneOpDim(storage, coef) {}

OperationLaplaceLinearBoundary::~OperationLaplaceLinearBoundary() {}

void OperationLaplaceLinearBoundary::up(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinearBoundary> s(func, *this->storage);
  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::down(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinearBoundary> s(func, *this->storage);
  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::downOpDim(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::upOpDim(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim) {
  UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
