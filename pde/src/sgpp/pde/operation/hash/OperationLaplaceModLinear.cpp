// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceModLinear.hpp>

#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp>
#include <sgpp/pde/basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplaceModLinear::OperationLaplaceModLinear(sgpp::base::GridStorage* storage)
    : UpDownOneOpDim(storage) {}

OperationLaplaceModLinear::~OperationLaplaceModLinear() {}

void OperationLaplaceModLinear::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                   size_t dim) {
  result.setAll(0.0);
  PhiPhiUpModLinear func(this->storage);
  sgpp::base::sweep<PhiPhiUpModLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) {
  result.setAll(0.0);
  PhiPhiDownModLinear func(this->storage);
  sgpp::base::sweep<PhiPhiDownModLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::downOpDim(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) {
  result.setAll(0.0);
  dPhidPhiDownModLinear func(this->storage);
  sgpp::base::sweep<dPhidPhiDownModLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::upOpDim(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) {
  result.setAll(0.0);
  dPhidPhiUpModLinear func(this->storage);
  sgpp::base::sweep<dPhidPhiUpModLinear> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
