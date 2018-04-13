// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretched.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp>
#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLTwoDotProductLinearStretched::OperationLTwoDotProductLinearStretched(
    sgpp::base::GridStorage* storage)
    : StdUpDown(storage) {}

OperationLTwoDotProductLinearStretched::~OperationLTwoDotProductLinearStretched() {}

void OperationLTwoDotProductLinearStretched::up(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiUpBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLTwoDotProductLinearStretched::down(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiDownBBLinearStretched func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinearStretched> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
