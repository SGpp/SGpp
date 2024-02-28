// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationSecondMomentLinearBoundary.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <cstdint>

namespace sgpp {
namespace base {

double OperationSecondMomentLinearBoundary::doQuadrature(DataVector& alpha, DataMatrix* bounds) {
  // handle bounds
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
        "OperationSecondMomentLinearBoundary::doQuadrature - bounds matrix has the wrong shape");
  }

  double res = 0;
  double tmpres = 1;
  double xlower = 0.0;
  double xupper = 0.0;
  double quad = 0.0;
  double firstMoment = 0.0;
  double secondMoment = 0.0;
  double value = 0.0;

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != storage.end(); iter++) {
    tmpres = 1.;

    for (size_t dim = 0; dim < storage.getDimension(); dim++) {
      auto index = iter->first->getIndex(dim);
      auto level = iter->first->getLevel(dim);
      double index_d = static_cast<double>(index);
      double level_d = static_cast<double>(level);

      // compute the second moment of basis function
      if (level == 0 && index == 0) {
        secondMoment = 0.25;
      } else if (level == 0 && index == 1) {
        secondMoment = 1. / 12.;
      } else {
        secondMoment = std::pow(8.0, -level_d) * (index_d * index_d + 1. / 6.);
      }

      if (bounds == nullptr) {
        // the second Moment is equal to the second moment of the basis function
        value = secondMoment;
      } else {
        // consider boundaries if given
        xlower = bounds->get(dim, 0);
        xupper = bounds->get(dim, 1);
        double width = (xupper - xlower);

        // compute the expected value of the basis function
        if (level == 0 && index == 0) {
          firstMoment = 1. / 3.;
        } else if (level == 0 && index == 1) {
          firstMoment = 1. / 6.;
        } else {
          firstMoment = index_d * std::pow(4.0, -level_d);
        }

        // compute the integral of the basis function
        if (level == 0) {
          quad = 0.5;
        } else {
          quad = std::pow(2.0, -level_d);
        }

        // compute the second moment of the transformed basis function
        value = width * width * secondMoment + 2 * width * xlower * firstMoment +
                xlower * xlower * quad;
      }

      tmpres *= value;
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
