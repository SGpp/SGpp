// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALMODWAVELETNAIVE_HPP
#define OPERATIONEVALMODWAVELETNAIVE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating modified wavelet linear combinations on Noboundary grids.
 */
class OperationEvalModWaveletNaive : public OperationEval {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   */
  explicit OperationEvalModWaveletNaive(GridStorage& storage) :
    storage(storage), pointInUnitCube(storage.getDimension()) {
  }

  /**
   * Destructor.
   */
  ~OperationEvalModWaveletNaive() override {
  }

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of the linear combination
   */
  double eval(const DataVector& alpha, const DataVector& point) override;

  /**
   * @param      alpha  coefficient matrix (each column is a coefficient vector)
   * @param      point  evaluation point
   * @param[out] value  values of linear combination
   */
  void eval(const DataMatrix& alpha, const DataVector& point,
            DataVector& value) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D wavelet basis
  SWaveletModifiedBase base;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALMODWAVELETNAIVE_HPP */
