// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

OperationLBLinear::OperationLBLinear(SGPP::base::GridStorage* storage)
    : SGPP::pde::StdUpDown(storage) {}

OperationLBLinear::~OperationLBLinear() {}

void OperationLBLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                           size_t dim) {
  // Dphi * phi
  DPhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLBLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                             size_t dim) {
  // Dphi * phi
  DPhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<DPhiPhiDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
