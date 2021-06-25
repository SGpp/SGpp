// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationQuadratureNakBsplineExtended.hpp"

#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

double OperationQuadratureNakBsplineExtended::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmpres = 0;

  for (size_t i = 0; i < storage.getSize(); i++) {
    tmpres = 1.;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *= base.getIntegral(storage.getPointLevel(i, d), storage.getPointIndex(i, d));
    }

    res += alpha[i] * tmpres;
  }

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  for (size_t d = 0; d < storage.getDimension(); d++) {
    res *= storage.getBoundingBox()->getIntervalWidth(d);
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
