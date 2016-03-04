// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationNaiveEvalPolyBoundary.hpp>

namespace sgpp {
namespace base {

double OperationNaiveEvalPolyBoundary::eval(const DataVector& alpha,
    const DataVector& point) {
  const size_t n = storage.getSize();
  const size_t dim = storage.getDimension();
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    const GridIndex& gp = *storage[i];
    double curValue = 1.0;

    for (size_t idim = 0; idim < dim; idim++) {
      const double val1d = base.evalSave(gp.getLevel(idim),
                                          gp.getIndex(idim), point[idim]);

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

}  // namespace base
}  // namespace sgpp
