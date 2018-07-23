// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationFirstMomentModLinear.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationFirstMomentModLinear::doQuadrature(const DataVector& alpha, DataMatrix* bounds) {
  // handle bounds
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
        "OperationFirstMomentModLinear::doQuadrature - bounds matrix has the wrong shape");
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
      double hInv = 1 << static_cast<int>(level);
      xlower = bounds == nullptr ? 0.0 : bounds->get(dim, 0);
      xupper = bounds == nullptr ? 1.0 : bounds->get(dim, 1);
      double width = xupper - xlower;
      if (level == 1)
        tmpres *= width * 0.5 + xlower;
      else if (index == 1)
        tmpres *= width * std::pow(4.0, 1 - level) / 3.0 + (2. / hInv) * xlower;
      else if (index == hInv - 1)
        tmpres *=
          width * (std::pow(2.0, 1 - 2*level) * (3. * hInv - 2) ) / 3. + (2. / hInv) * xlower;
      else
        tmpres *=
          width * index * std::pow(4.0, -level) + xlower * std::pow(2.0, -level);
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
