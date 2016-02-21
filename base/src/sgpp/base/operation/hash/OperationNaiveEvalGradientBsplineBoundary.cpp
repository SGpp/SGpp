// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientBsplineBoundary.hpp>

namespace SGPP {
namespace base {

float_t OperationNaiveEvalGradientBsplineBoundary::evalGradient(
  const DataVector& alpha, const DataVector& point, DataVector& gradient) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  float_t result = 0.0;

  gradient.resize(storage.getDimension());
  gradient.setAll(0.0);

  DataVector curGradient(d);

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    float_t curValue = 1.0;
    curGradient.setAll(alpha[i]);

    for (size_t t = 0; t < d; t++) {
      const float_t val1d = base.eval(gp.getLevel(t), gp.getIndex(t), point[t]);
      const float_t dx1d = base.evalDx(gp.getLevel(t), gp.getIndex(t),
                                       point[t]);

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

}  // namespace base
}  // namespace SGPP
