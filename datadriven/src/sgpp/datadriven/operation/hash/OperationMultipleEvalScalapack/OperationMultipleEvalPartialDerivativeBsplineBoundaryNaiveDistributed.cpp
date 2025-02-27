/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalPartialDerivativeBsplineBoundaryNaiveDistributed.cpp
 *
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalPartialDerivativeBsplineBoundaryNaiveDistributed.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

void OperationMultipleEvalPartialDerivativeBsplineBoundaryNaiveDistributed::multDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  if (!result.isProcessMapped()) {
    return;
  }

  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  const double innerDerivative = 1.0 / storage.getBoundingBox()->getIntervalWidth(derivDim);

  for (size_t localJ = 0; localJ < result.getLocalRows(); ++localJ) {
    size_t j = result.localToGlobalRowIndex(localJ);
    for (size_t i = 0; i < n; ++i) {
      const base::GridPoint& gp = storage[i];
      double curPartialDerivative = 1.0;

      const double dx1d =
          base.evalDx(gp.getLevel(derivDim), gp.getIndex(derivDim), pointsInUnitCube(j, derivDim)) *
          innerDerivative;

      for (size_t t = 0; t < d; ++t) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if ((val1d == 0.0) && (t != derivDim)) {
          curPartialDerivative = 0.0;
          break;
        }

        curPartialDerivative *= ((t == derivDim) ? dx1d : val1d);
      }

      double newResult = result.get(j) + alpha[i] * curPartialDerivative;

      result.set(j, newResult);
    }
  }
}

void OperationMultipleEvalPartialDerivativeBsplineBoundaryNaiveDistributed::
    multTransposeDistributed(DataVector& alpha, DataVectorDistributed& result) {
  if (!result.isProcessMapped()) {
    return;
  }

  const size_t d = storage.getDimension();
  const size_t m = dataset.getNrows();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  const double innerDerivative = 1.0 / storage.getBoundingBox()->getIntervalWidth(derivDim);

  for (size_t iLocal = 0; iLocal < result.getLocalRows(); ++iLocal) {
    size_t i = result.localToGlobalRowIndex(iLocal);
    const base::GridPoint& gp = storage[i];

    for (size_t j = 0; j < m; ++j) {
      double curPartialDerivative = 1.0;

      const double dx1d =
          base.evalDx(gp.getLevel(derivDim), gp.getIndex(derivDim), pointsInUnitCube(j, derivDim)) *
          innerDerivative;

      for (size_t t = 0; t < d; ++t) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if ((val1d == 0.0) && (t != derivDim)) {
          curPartialDerivative = 0.0;
          break;
        }

        curPartialDerivative *= ((t == derivDim) ? dx1d : val1d);
      }

      double newResult = result.get(i) + alpha[j] * curPartialDerivative;

      result.set(i, newResult);
    }
  }
}

double OperationMultipleEvalPartialDerivativeBsplineBoundaryNaiveDistributed::getDuration() {
  return 0.0;
}

}  // namespace datadriven
}  // namespace sgpp
