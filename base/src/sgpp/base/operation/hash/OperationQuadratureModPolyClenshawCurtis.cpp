// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureModPolyClenshawCurtis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationQuadratureModPolyClenshawCurtis::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmpres = 0;

  for (size_t i = 0; i < alpha.getSize(); i++) {
    GridPoint& gp = storage.getPoint(i);
    tmpres = 1.;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *= base.getIntegral(gp.getLevel(d), gp.getIndex(d));
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
