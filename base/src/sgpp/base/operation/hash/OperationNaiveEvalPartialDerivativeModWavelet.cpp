// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeModWavelet.hpp>

namespace sgpp {
namespace base {

double OperationNaiveEvalPartialDerivativeModWavelet::evalPartialDerivative(
    const DataVector& alpha,
    const DataVector& point,
    size_t derivDim) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    double curValue = 1.0;

    for (size_t t = 0; t < d; t++) {
      const double val1d = ((t == derivDim) ?
                             base.evalDx(gp.getLevel(t), gp.getIndex(t), point[t]) :
                             base.eval(gp.getLevel(t), gp.getIndex(t), point[t]));

      if (val1d == 0.0) {
        curValue = 0.0;
        break;
      }

      curValue *= val1d;
    }

    result += alpha[i] * curValue;
  }

  return result;
}

void OperationNaiveEvalPartialDerivativeModWavelet::evalPartialDerivative(
    const DataMatrix& alpha,
    const DataVector& point,
    size_t derivDim,
    DataVector& partialDerivative) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = alpha.getNcols();

  partialDerivative.resize(m);
  partialDerivative.setAll(0.0);

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    double curValue = 1.0;

    for (size_t t = 0; t < d; t++) {
      const double val1d = ((t == derivDim) ?
                             base.evalDx(gp.getLevel(t), gp.getIndex(t), point[t]) :
                             base.eval(gp.getLevel(t), gp.getIndex(t), point[t]));

      if (val1d == 0.0) {
        curValue = 0.0;
        break;
      }

      curValue *= val1d;
    }

    for (size_t j = 0; j < m; j++) {
      partialDerivative[j] += alpha(i, j) * curValue;
    }
  }
}

}  // namespace base
}  // namespace sgpp
