// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationSecondMomentModPolyClenshawCurtis.hpp>
#include <sgpp/base/grid/type/ModPolyClenshawCurtisGrid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationSecondMomentModPolyClenshawCurtis::doQuadrature(DataVector& alpha,
                                                                DataMatrix* bounds) {
  // handle bounds
  GridStorage& storage = grid->getStorage();
  size_t numDims = storage.getDimension();

  // check if the boundaries are given in the right shape
  if (bounds != nullptr && (bounds->getNcols() != 2 || bounds->getNrows() != numDims)) {
    throw application_exception(
          "OperationSecondMomentPoly::doQuadrature - bounds matrix has the wrong shape");
  }

  double res = 0;
  double tmpres = 1;
  base::index_t index;
  double indexDbl;
  base::level_t level;
  double xlower = 0.0;
  double xupper = 0.0;
  const size_t quadOrder =  static_cast<size_t>(
    ceil(static_cast<double>(
         dynamic_cast<sgpp::base::ModPolyClenshawCurtisGrid*>(grid)->getDegree()) / 2.))
    + 2;
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
      double left = (index == 0) ? 0.0 : clenshawCurtisTable.getPoint(level, index - 1);
      double right = (indexDbl == hInv) ? 1.0 : clenshawCurtisTable.getPoint(level, index + 1);
      double scaling = right - left;

      double gaussQuadSumSecondMoment = 0.;
      double gaussQuadSumFirstMoment = 0.;
      for (size_t c = 0; c < quadOrder; c++) {
        const double x = left + scaling * coordinates[c];
        gaussQuadSumSecondMoment += weights[c] * x * x * basis.eval(level, index, x);
        gaussQuadSumFirstMoment += weights[c] * x * basis.eval(level, index, x);
      }

      gaussQuadSumSecondMoment *= scaling;
      gaussQuadSumFirstMoment *= scaling;
      tmpres *=
        width * width * gaussQuadSumSecondMoment
        + 2 * width * xlower * gaussQuadSumFirstMoment
        + xlower * xlower * basis.getIntegral(level, index);
    }

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
