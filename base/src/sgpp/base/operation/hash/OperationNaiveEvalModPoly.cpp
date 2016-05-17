// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationNaiveEvalModPoly.hpp>

namespace sgpp {
namespace base {

double OperationNaiveEvalModPoly::eval(const DataVector& alpha, const DataVector& point) {
  const size_t n = storage.getSize();
  const size_t dim = storage.getDimension();
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    double curValue = 1.0;

    for (size_t idim = 0; idim < dim; idim++) {
      const double val1d = base.evalSave(gp.getLevel(idim), gp.getIndex(idim), point[idim]);

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

void OperationNaiveEvalModPoly::eval(const DataMatrix& alpha, const DataVector& point,
                                     DataVector& value) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = alpha.getNcols();

  value.resize(m);
  value.setAll(0.0);

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    double curValue = 1.0;

    for (size_t t = 0; t < d; t++) {
      const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), point[t]);

      if (val1d == 0.0) {
        curValue = 0.0;
        break;
      }

      curValue *= val1d;
    }

    for (size_t j = 0; j < m; j++) {
      value[j] += alpha(i, j) * curValue;
    }
  }
}

}  // namespace base
}  // namespace sgpp
