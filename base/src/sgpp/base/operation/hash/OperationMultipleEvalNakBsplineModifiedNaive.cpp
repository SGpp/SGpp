// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationMultipleEvalNakBsplineModifiedNaive.hpp"

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationMultipleEvalNakBsplineModifiedNaive::mult(DataVector& alpha,
                                                             DataVector& result) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = dataset.getNrows();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      const GridPoint& gp = storage[i];
      double curValue = 1.0;

      for (size_t t = 0; t < d; t++) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if (val1d == 0.0) {
          curValue = 0.0;
          break;
        }

        curValue *= val1d;
      }

      result[j] += alpha[i] * curValue;
    }
  }
}

void OperationMultipleEvalNakBsplineModifiedNaive::multTranspose(DataVector& alpha,
                                                                      DataVector& result) {
  const size_t n = storage.getSize();
  const size_t d = storage.getDimension();
  const size_t m = dataset.getNrows();

  result.setAll(0.0);

  pointsInUnitCube = dataset;
  storage.getBoundingBox()->transformPointsToUnitCube(pointsInUnitCube);

  for (size_t i = 0; i < n; i++) {
    const GridPoint& gp = storage[i];

    for (size_t j = 0; j < m; j++) {
      double curValue = 1.0;

      for (size_t t = 0; t < d; t++) {
        const double val1d = base.eval(gp.getLevel(t), gp.getIndex(t), pointsInUnitCube(j, t));

        if (val1d == 0.0) {
          curValue = 0.0;
          break;
        }

        curValue *= val1d;
      }

      result[i] += alpha[j] * curValue;
    }
  }
}

double OperationMultipleEvalNakBsplineModifiedNaive::getDuration() { return 0.0; }

}  // namespace base
}  // namespace sgpp
