// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationSecondMomentLinear.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationSecondMomentLinear::doQuadrature(DataVector& alpha, DataMatrix* bounds) {
  // handle bounds
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
        "OperationSecondMomentLinear::doQuadrature - bounds matrix has the wrong shape");
  }

  double res = 0;
  double tmpres = 1;
  double index;
  double level;
  double xlower = 0.0;
  double xupper = 0.0;

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != storage.end(); iter++) {
    tmpres = 1.;

    for (size_t dim = 0; dim < storage.getDimension(); dim++) {
      index = static_cast<double>(iter->first->getIndex(dim));
      level = static_cast<double>(iter->first->getLevel(dim));

      if (bounds == nullptr) {
        // evaluate the second moment in the unit hypercube
        tmpres *= pow(8.0, -level) * (index * index + 1. / 6.);
      } else {
        // consider boundaries if given
        xlower = bounds->get(dim, 0);
        xupper = bounds->get(dim, 1);
        double width = (xupper - xlower);
        tmpres *= width * width * (index * index + 1. / 6.) * std::pow(8, -level) +
                  2 * width * xlower * index * std::pow(4.0, -level) +
                  xlower * xlower * std::pow(2.0, -level);
      }
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
