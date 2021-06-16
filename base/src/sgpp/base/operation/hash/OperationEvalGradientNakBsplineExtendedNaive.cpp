// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationEvalGradientNakBsplineExtendedNaive.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationEvalGradientNakBsplineExtendedNaive::evalGradient(const DataVector& alpha,
                                                                  const DataVector& point,
                                                                  DataVector& gradient) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  double result = 0.0;

  pointInUnitCube = point;
  storage.getBoundingBox()->transformPointToUnitCube(pointInUnitCube);

  for (size_t t = 0; t < d; t++) {
    innerDerivative[t] = 1.0 / storage.getBoundingBox()->getIntervalWidth(t);
  }

  gradient.resize(storage.getDimension());
  gradient.setAll(0.0);

  DataVector curGradient(d);

  for (size_t i = 0; i < n; i++) {
    const GridPoint& gp = storage[i];
    double curValue = 1.0;
    curGradient.setAll(alpha[i]);

    for (size_t t = 0; t < d; t++) {
      const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]);
      const double dx1d =
          base.evalDx(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]) * innerDerivative[t];

      curValue *= val1d;

      for (size_t t2 = 0; t2 < d; t2++) {
        if (t2 == t) {
          curGradient[t2] *= dx1d;
        } else {
          curGradient[t2] *= val1d;
        }
      }
    }

    result += alpha[i] * curValue;
    gradient.add(curGradient);
  }

  return result;
}

void OperationEvalGradientNakBsplineExtendedNaive::evalGradient(const DataMatrix& alpha,
                                                                const DataVector& point,
                                                                DataVector& value,
                                                                DataMatrix& gradient) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = alpha.getNcols();

  pointInUnitCube = point;
  storage.getBoundingBox()->transformPointToUnitCube(pointInUnitCube);

  for (size_t t = 0; t < d; t++) {
    innerDerivative[t] = 1.0 / storage.getBoundingBox()->getIntervalWidth(t);
  }

  value.resize(m);
  value.setAll(0.0);

  gradient.resize(m, d);
  gradient.setAll(0.0);

  DataVector curGradient(d);

  for (size_t i = 0; i < n; i++) {
    const GridPoint& gp = storage[i];
    double curValue = 1.0;
    curGradient.setAll(1.0);

    for (size_t t = 0; t < d; t++) {
      const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]);
      const double dx1d =
          base.evalDx(gp.getLevel(t), gp.getIndex(t), pointInUnitCube[t]) * innerDerivative[t];

      curValue *= val1d;

      for (size_t t2 = 0; t2 < d; t2++) {
        if (t2 == t) {
          curGradient[t2] *= dx1d;
        } else {
          curGradient[t2] *= val1d;
        }
      }
    }

    for (size_t j = 0; j < m; j++) {
      value[j] += alpha(i, j) * curValue;

      for (size_t t = 0; t < d; t++) {
        gradient(j, t) += alpha(i, j) * curGradient[t];
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp
