// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationWeightedQuadratureNakBsplineModified.hpp"

namespace sgpp {
namespace base {

double OperationWeightedQuadratureNakBsplineModified::doWeightedQuadrature(
    DataVector& alpha, std::shared_ptr<sgpp::base::Distribution> pdf, size_t quadOrder) {
  double mean = 0;
  for (size_t i = 0; i < storage.getSize(); i++) {
    double tmpres = 1;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      tmpres *=
          base.getMean(storage.getPointLevel(i, d), storage.getPointIndex(i, d), pdf, quadOrder);
    }
    mean += alpha[i] * tmpres;
  }
  return mean;
}

}  // namespace base
}  // namespace sgpp
