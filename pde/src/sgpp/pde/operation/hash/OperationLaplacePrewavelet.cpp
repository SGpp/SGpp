// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceDownGradientPrewavelet.hpp>
#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpGradientPrewavelet.hpp>
#include <sgpp/pde/basis/prewavelet/algorithm_sweep/LaplaceUpPrewavelet.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePrewavelet.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLaplacePrewavelet::OperationLaplacePrewavelet(SGPP::base::GridStorage* storage, SGPP::base::GridStorage* shadowstorage) :
      UpDownOneOpDimWithShadow(storage, shadowstorage) {

    }

    OperationLaplacePrewavelet::~OperationLaplacePrewavelet() {
    }

    void OperationLaplacePrewavelet::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      LaplaceUpPrewavelet func(this->storage);
      SGPP::base::sweep<LaplaceUpPrewavelet> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplacePrewavelet::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
    }

    void OperationLaplacePrewavelet::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      LaplaceDownGradientPrewavelet func(this->storage);
      SGPP::base::sweep<LaplaceDownGradientPrewavelet> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplacePrewavelet::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      LaplaceUpGradientPrewavelet func(this->storage);
      SGPP::base::sweep<LaplaceUpGradientPrewavelet> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

  }
}