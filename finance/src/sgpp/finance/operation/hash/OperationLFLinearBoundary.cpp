// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLFLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationLFLinearBoundary::OperationLFLinearBoundary(sgpp::base::GridStorage* storage)
    : sgpp::pde::StdUpDown(storage) {}

OperationLFLinearBoundary::~OperationLFLinearBoundary() {}

void OperationLFLinearBoundary::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                   size_t dim) {
  // X * dphi * phi
  XdPhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLFLinearBoundary::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {
  // X * dphi * phi
  XdPhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
