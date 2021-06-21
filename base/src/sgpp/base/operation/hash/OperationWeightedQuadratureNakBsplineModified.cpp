// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationWeightedQuadratureNakBsplineModified.hpp"

namespace sgpp {
namespace base {

double OperationWeightedQuadratureNakBsplineModified::doWeightedQuadrature(
    DataVector& alpha, sgpp::base::DistributionsVector pdfs) {
  double mean = 0;
  for (size_t i = 0; i < storage.getSize(); i++) {
    double tmpres = 1;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *= base.getMean(storage.getPointLevel(i, d), storage.getPointIndex(i, d), pdfs.get(d),
                             quadCoordinates, quadWeights);
    }
    mean += alpha[i] * tmpres;
  }
  return mean;
}

}  // namespace base
}  // namespace sgpp
