// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLDLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

OperationLDLinearBoundary::OperationLDLinearBoundary(SGPP::base::GridStorage* storage)
    : SGPP::pde::StdUpDown(storage) {}

OperationLDLinearBoundary::~OperationLDLinearBoundary() {}

void OperationLDLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                   size_t dim) {
  // X * phi * phi
  XPhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLDLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                     size_t dim) {
  // X * phi * phi
  XPhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
