// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadraturePoly.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationQuadraturePoly::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmpres = 0;
  GridIndex* gp;

  for (size_t i = 0; i < alpha.getSize(); i++) {
    gp = storage.getGridIndex(i);
    tmpres = 1.;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *= base.getIntegral(gp->getLevel(d), gp->getIndex(d));
    }

    res += alpha[i] * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
