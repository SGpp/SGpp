// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceDownGradientPrewavelet.hpp>
#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpGradientPrewavelet.hpp>
#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpPrewavelet.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePrewavelet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

OperationLaplacePrewavelet::OperationLaplacePrewavelet(sgpp::base::GridStorage* storage,
                                                       sgpp::base::GridStorage* shadowstorage)
    : UpDownOneOpDimWithShadow(storage, shadowstorage) {}

OperationLaplacePrewavelet::~OperationLaplacePrewavelet() {}

void OperationLaplacePrewavelet::up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                    size_t dim) {
  LaplaceUpPrewavelet func(this->storage);
  sgpp::base::sweep<LaplaceUpPrewavelet> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplacePrewavelet::down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) {}

void OperationLaplacePrewavelet::downOpDim(sgpp::base::DataVector& alpha,
                                           sgpp::base::DataVector& result, size_t dim) {
  LaplaceDownGradientPrewavelet func(this->storage);
  sgpp::base::sweep<LaplaceDownGradientPrewavelet> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}

void OperationLaplacePrewavelet::upOpDim(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) {
  LaplaceUpGradientPrewavelet func(this->storage);
  sgpp::base::sweep<LaplaceUpGradientPrewavelet> s(func, *this->storage);
  s.sweep1D(alpha, result, dim);
}
}  // namespace pde
}  // namespace sgpp
