// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLTwoDotProductLinearBoundary::OperationLTwoDotProductLinearBoundary(
    sgpp::base::GridStorage* storage)
    : StdUpDown(storage) {}

OperationLTwoDotProductLinearBoundary::~OperationLTwoDotProductLinearBoundary() {}

void OperationLTwoDotProductLinearBoundary::up(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLTwoDotProductLinearBoundary::down(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
