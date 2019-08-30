// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeModNakBsplineNaive.hpp>

namespace sgpp {
namespace base {

double OperationEvalPartialDerivativeModNakBsplineNaive::evalPartialDerivative(
    const DataVector& alpha,
    const DataVector& point,
    size_t derivDim,
    double& partialDerivative) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  double result = 0.0;

  partialDerivative = 0.0;
  pointInUnitCube = point;
  storage.getBoundingBox()->transformPointToUnitCube(pointInUnitCube);

  const double innerDerivative = 1.0 / storage.getBoundingBox()->getIntervalWidth(derivDim);

  for (size_t i = 0; i < n; i++) {
    const GridPoint& gp = storage[i];
    double curValue = 1.0;
    double curPartialDerivative = 1.0;

    const double dx1d = baseDeriv1.eval(gp.getLevel(derivDim), gp.getIndex(derivDim),
                                    pointInUnitCube[derivDim]) * innerDerivative;

    for (size_t t = 0; t < d; t++) {
      const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]);

      if ((val1d == 0.0) && (t != derivDim)) {
        curValue = 0.0;
        curPartialDerivative = 0.0;
        break;
      }

      curValue *= val1d;
      curPartialDerivative *= ((t == derivDim) ? dx1d : val1d);
    }

    result += alpha[i] * curValue;
    partialDerivative += alpha[i] * curPartialDerivative;
  }

  return result;
}

void OperationEvalPartialDerivativeModNakBsplineNaive::evalPartialDerivative(
    const DataMatrix& alpha,
    const DataVector& point,
    size_t derivDim,
    DataVector& value,
    DataVector& partialDerivative) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = alpha.getNcols();

  pointInUnitCube = point;
  storage.getBoundingBox()->transformPointToUnitCube(pointInUnitCube);

  const double innerDerivative = 1.0 / storage.getBoundingBox()->getIntervalWidth(derivDim);

  value.resize(m);
  value.setAll(0.0);
  partialDerivative.resize(m);
  partialDerivative.setAll(0.0);

  for (size_t i = 0; i < n; i++) {
    const GridPoint& gp = storage[i];
    double curValue = 1.0;
    double curPartialDerivative = 1.0;

    const double dx1d = baseDeriv1.eval(gp.getLevel(derivDim), gp.getIndex(derivDim),
                                    pointInUnitCube[derivDim]) * innerDerivative;

    for (size_t t = 0; t < d; t++) {
      const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]);

      if ((val1d == 0.0) && (t != derivDim)) {
        curValue = 0.0;
        curPartialDerivative = 0.0;
        break;
      }

      curValue *= val1d;
      curPartialDerivative *= ((t == derivDim) ? dx1d : val1d);
    }

    for (size_t j = 0; j < m; j++) {
      value[j] += alpha(i, j) * curValue;
      partialDerivative[j] += alpha(i, j) * curPartialDerivative;
    }
  }
}

}  // namespace base
}  // namespace sgpp
