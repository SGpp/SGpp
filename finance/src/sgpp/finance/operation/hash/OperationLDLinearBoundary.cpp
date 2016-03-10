// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLDLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationLDLinearBoundary::OperationLDLinearBoundary(sgpp::base::GridStorage* storage)
    : sgpp::pde::StdUpDown(storage) {}

OperationLDLinearBoundary::~OperationLDLinearBoundary() {}

void OperationLDLinearBoundary::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                   size_t dim) {
  // X * phi * phi
  XPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLDLinearBoundary::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {
  // X * phi * phi
  XPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
