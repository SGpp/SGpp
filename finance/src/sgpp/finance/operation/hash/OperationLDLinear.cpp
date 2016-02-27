// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLDLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

OperationLDLinear::OperationLDLinear(sgpp::base::GridStorage* storage)
    : sgpp::pde::StdUpDown(storage) {}

OperationLDLinear::~OperationLDLinear() {}

void OperationLDLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                           size_t dim) {
  // X * phi * phi
  XPhiPhiUpBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLDLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                             size_t dim) {
  // X * phi * phi
  XPhiPhiDownBBLinear func(this->storage);
  sgpp::base::sweep<XPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
