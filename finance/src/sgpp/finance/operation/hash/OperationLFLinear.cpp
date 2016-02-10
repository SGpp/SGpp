// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLFLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

OperationLFLinear::OperationLFLinear(SGPP::base::GridStorage* storage) :
  SGPP::pde::StdUpDown(storage) {
}

OperationLFLinear::~OperationLFLinear() {
}

void OperationLFLinear::up(SGPP::base::DataVector& alpha,
                           SGPP::base::DataVector& result, size_t dim) {
  // X * dphi * phi
  XdPhiPhiUpBBLinear func(this->storage);
  SGPP::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLFLinear::down(SGPP::base::DataVector& alpha,
                             SGPP::base::DataVector& result, size_t dim) {
  // X * dphi * phi
  XdPhiPhiDownBBLinear func(this->storage);
  SGPP::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

  s.sweep1D(alpha, result, dim);
}

}
}