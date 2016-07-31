// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureLinearBoundary.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

double OperationQuadratureLinearBoundary::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmp;
  int nr_boundaries = 0;  // nr. of boundaries a point touches
  int cur_ind, cur_lev;
  GridStorage::point_type index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    tmp = pow(2.0, -static_cast<double>(iter->first->getLevelSum())) * alpha.get(
            iter->second);

    if (!iter->first->isInnerPoint()) {
      nr_boundaries = 0;

      for (size_t d = 0; d < iter->first->getDimension(); d++) {
        cur_ind = iter->first->getIndex(d);
        cur_lev = iter->first->getLevel(d);

        if ((cur_ind == 0) || (pow(2.0, cur_lev) == cur_ind))
          nr_boundaries++;
      }

      tmp *= (pow(2.0, -nr_boundaries));
    }

    res += tmp;
  }

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  for (size_t d = 0; d < storage.getDimension(); d++) {
    res *= storage.getBoundingBox()->getIntervalWidth(d);
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
