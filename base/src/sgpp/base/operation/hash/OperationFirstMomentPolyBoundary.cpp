// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationFirstMomentPolyBoundary.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationFirstMomentPolyBoundary::doQuadrature(const DataVector& alpha, DataMatrix* bounds) {
  // handle bounds
  GridStorage& storage = grid->getStorage();
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
        "OperationFirstMomentPolyBoundary::doQuadrature - bounds matrix has the wrong shape");
  }

  double res = 0;
  double tmpres = 1;
  base::index_t index;
  double indexDbl;
  base::level_t level;
  double xlower = 0.0;
  double xupper = 0.0;
  const size_t quadOrder =  static_cast<size_t>(
    ceil(static_cast<double>(dynamic_cast<sgpp::base::PolyBoundaryGrid*>(grid)->getDegree()) / 2.))
    + 1;
  base::SBasis& basis = const_cast<base::SBasis&>(grid->getBasis());
  base::DataVector coordinates;
  base::DataVector weights;
  base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != storage.end(); iter++) {
    tmpres = 1.;

    for (size_t dim = 0; dim < storage.getDimension(); dim++) {
      index = iter->first->getIndex(dim);
      indexDbl = static_cast<double>(index);
      level = iter->first->getLevel(dim);
      double hInv = static_cast<double>(1 << level);
      xlower = bounds == nullptr ? 0.0 : bounds->get(dim, 0);
      xupper = bounds == nullptr ? 1.0 : bounds->get(dim, 1);
      double width = xupper - xlower;
      double scaling = (index != 0 && indexDbl != hInv) ? 2./hInv : 1./hInv;
      double left = (index != 0) ? (indexDbl - 1) * (1./hInv) : 0.0;

      double gaussQuadSum = 0.;
      for (size_t c = 0; c < quadOrder; c++) {
        const double x = left + scaling * coordinates[c];
        gaussQuadSum += weights[c] * x * basis.eval(level, index, x);
      }

      tmpres *=
        width * scaling * gaussQuadSum + xlower * basis.getIntegral(level, index);
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
