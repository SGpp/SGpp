// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceEnhancedLinearBoundary::OperationLaplaceEnhancedLinearBoundary(
    sgpp::base::GridStorage* storage)
    : UpDownOneOpDimEnhanced(storage) {}

OperationLaplaceEnhancedLinearBoundary::OperationLaplaceEnhancedLinearBoundary(
    sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef)
    : UpDownOneOpDimEnhanced(storage, coef) {}

OperationLaplaceEnhancedLinearBoundary::~OperationLaplaceEnhancedLinearBoundary() {}

void OperationLaplaceEnhancedLinearBoundary::up(sgpp::base::DataMatrix& alpha,
                                                sgpp::base::DataMatrix& result, size_t dim) {
  LaplaceEnhancedUpBBLinearBoundary func(this->storage);
  sgpp::base::sweep<LaplaceEnhancedUpBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceEnhancedLinearBoundary::down(sgpp::base::DataMatrix& alpha,
                                                  sgpp::base::DataMatrix& result, size_t dim) {
  LaplaceEnhancedDownBBLinearBoundary func(this->storage);
  sgpp::base::sweep<LaplaceEnhancedDownBBLinearBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
