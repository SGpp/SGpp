// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

OperationLTwoDotProductLinearBoundary::OperationLTwoDotProductLinearBoundary(
    SGPP::base::GridStorage* storage)
    : StdUpDown(storage) {}

OperationLTwoDotProductLinearBoundary::~OperationLTwoDotProductLinearBoundary() {}

void OperationLTwoDotProductLinearBoundary::up(SGPP::base::DataVector& alpha,
                                               SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiUpBBLinearBoundary func(this->storage);
  SGPP::base::sweep<PhiPhiUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLTwoDotProductLinearBoundary::down(SGPP::base::DataVector& alpha,
                                                 SGPP::base::DataVector& result, size_t dim) {
  // phi * phi
  PhiPhiDownBBLinearBoundary func(this->storage);
  SGPP::base::sweep<PhiPhiDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace pde
}  // namespace SGPP
