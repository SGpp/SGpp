// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationFirstMomentLinearBoundary.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <cstdint>

namespace sgpp {
namespace base {

double OperationFirstMomentLinearBoundary::doQuadrature(const DataVector& alpha,
                                                        DataMatrix* bounds) {
  // handle bounds
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
        "OperationFirstMomentLinearBoundary::doQuadrature - bounds matrix has the wrong shape");
  }

  double res = 0;
  double tmpres = 1;
  double xlower = 0.0;
  double xupper = 0.0;
  double value = 0.0;
  double mean = 0.0;
  double quad = 0.0;

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != storage.end(); iter++) {
    tmpres = 1.;

    for (size_t dim = 0; dim < storage.getDimension(); dim++) {
      auto index = iter->first->getIndex(dim);
      auto level = iter->first->getLevel(dim);
      double index_d = static_cast<double>(index);
      double level_d = static_cast<double>(level);

      // compute the expected value of the basis function
      if (level == 0) {
        if (index == 0) {
          mean = 1. / 3.;
        } else {
          mean = 1. / 6.;
        }
      } else {
        // evaluate the first moment in the unit hypercube
        mean = index_d * std::pow(4.0, -level_d);
      }

      // compute the integral of the basis function if necessary
      if (bounds != nullptr) {
        // consider boundaries if given
        xlower = bounds->get(dim, 0);
        xupper = bounds->get(dim, 1);
        if (level == 0) {
          quad = 0.5;
        } else {
          quad = std::pow(2.0, -level_d);
        }

        // compute the first moment for the current basis function
        value = (xupper - xlower) * mean + xlower * quad;
      } else {
        // for the unit interval the first moment is equal to the mean
        value = mean;
      }

      tmpres *= value;
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
