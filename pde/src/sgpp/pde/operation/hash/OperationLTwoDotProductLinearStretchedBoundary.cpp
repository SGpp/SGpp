// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLTwoDotProductLinearStretchedBoundary::OperationLTwoDotProductLinearStretchedBoundary(
    sgpp::base::GridStorage* storage)
    : StdUpDown(storage) {}

OperationLTwoDotProductLinearStretchedBoundary::~OperationLTwoDotProductLinearStretchedBoundary() {}

void OperationLTwoDotProductLinearStretchedBoundary::up(sgpp::base::DataVector& alpha,
                                                        sgpp::base::DataVector& result,
                                                        size_t dim) {
  // phi * phi
  PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLTwoDotProductLinearStretchedBoundary::down(sgpp::base::DataVector& alpha,
                                                          sgpp::base::DataVector& result,
                                                          size_t dim) {
  // phi * phi
  PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
