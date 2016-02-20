// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadraturePolyBoundary.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

float_t OperationQuadraturePolyBoundary::doQuadrature(DataVector& alpha) {
  float_t res = 0;
  float_t tmpres = 0;
  GridIndex* gp;

  for (size_t i = 0; i < alpha.getSize(); i++) {
    gp = storage.get(i);
    tmpres = 1.;

    for (size_t d = 0; d < storage.dim(); d++) {
      tmpres *= base.getIntegral(gp->getLevel(d), gp->getIndex(d));
    }

    res += alpha[i] * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace SGPP
