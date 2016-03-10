// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationLBLinear::OperationLBLinear(sgpp::base::GridStorage* storage)
    : sgpp::pde::StdUpDown(storage) {}

OperationLBLinear::~OperationLBLinear() {}

void OperationLBLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                           size_t dim) {
  // Dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLBLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                             size_t dim) {
  // Dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
