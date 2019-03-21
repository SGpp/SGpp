// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationWeightedQuadratureNakBsplineModifiedCombigrid.hpp"

#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace combigrid {

double OperationWeightedQuadratureNakBsplineModifiedCombigrid::doQuadrature(
    sgpp::base::DataVector& alpha) {
  double res = 0;
  double tmpres = 0;

  size_t numAdditionalPoints = 0;
  size_t incrementQuadraturePoints = 10;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    sgpp::base::GridPoint& gp = storage.getPoint(i);
    tmpres = 1.;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *= base.getWeightedIntegral(
          gp.getLevel(d), gp.getIndex(d), weightFunctionsCollection[d], bounds[2 * d],
          bounds[2 * d + 1], numAdditionalPoints, incrementQuadraturePoints);
    }

    res += alpha[i] * tmpres;
  }

  // This is done in other OperationQuadrature routines. Necessary here?
  // multiply with determinant of "unit cube -> BoundingBox" transformation
  //  for (size_t d = 0; d < storage.getDimension(); d++) {
  //    res *= storage.getBoundingBox()->getIntervalWidth(d);
  //  }

  return res;
}

}  // namespace combigrid
}  // namespace sgpp
