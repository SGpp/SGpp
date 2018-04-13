// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(sgpp::base::GridStorage* storage)
    : UpDownOneOpDimEnhanced(storage) {}

OperationLaplaceEnhancedLinear::OperationLaplaceEnhancedLinear(sgpp::base::GridStorage* storage,
                                                               sgpp::base::DataVector& coef)
    : UpDownOneOpDimEnhanced(storage, coef) {}

OperationLaplaceEnhancedLinear::~OperationLaplaceEnhancedLinear() {}

void OperationLaplaceEnhancedLinear::up(sgpp::base::DataMatrix& alpha,
                                        sgpp::base::DataMatrix& result, size_t dim) {
  LaplaceEnhancedUpBBLinear func(this->storage);
  sgpp::base::sweep<LaplaceEnhancedUpBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceEnhancedLinear::down(sgpp::base::DataMatrix& alpha,
                                          sgpp::base::DataMatrix& result, size_t dim) {
  LaplaceEnhancedDownBBLinear func(this->storage);
  sgpp::base::sweep<LaplaceEnhancedDownBBLinear> s(func, *this->storage);

  s.sweep1D(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
