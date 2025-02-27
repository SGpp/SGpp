/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalModBsplineNaiveDistributed.cpp
 *
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalModBsplineNaiveDistributed.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

void OperationMultipleEvalModBsplineNaiveDistributed::multDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  if (!result.isProcessMapped()) {
    return;
  }

  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  for (size_t localJ = 0; localJ < result.getLocalRows(); ++localJ) {
    size_t j = result.localToGlobalRowIndex(localJ);
    for (size_t i = 0; i < n; ++i) {
      const base::GridPoint& gp = storage[i];
      double curValue = 1.0;

      for (size_t t = 0; t < d; ++t) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if (val1d == 0.0) {
          curValue = 0.0;
          break;
        }

        curValue *= val1d;
      }

      double newResult = result.get(j) + alpha[i] * curValue;

      result.set(j, newResult);
    }
  }
}

void OperationMultipleEvalModBsplineNaiveDistributed::multTransposeDistributed(
    DataVector& alpha, DataVectorDistributed& result) {
  if (!result.isProcessMapped()) {
    return;
  }

  const size_t d = storage.getDimension();
  const size_t m = dataset.getNrows();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  for (size_t iLocal = 0; iLocal < result.getLocalRows(); ++iLocal) {
    size_t i = result.localToGlobalRowIndex(iLocal);
    const base::GridPoint& gp = storage[i];

    for (size_t j = 0; j < m; ++j) {
      double curValue = 1.0;

      for (size_t t = 0; t < d; ++t) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if (val1d == 0.0) {
          curValue = 0.0;
          break;
        }

        curValue *= val1d;
      }

      double newResult = result.get(i) + alpha[j] * curValue;

      result.set(i, newResult);
    }
  }
}

double OperationMultipleEvalModBsplineNaiveDistributed::getDuration() { return 0.0; }

}  // namespace datadriven
}  // namespace sgpp
