// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureModLinear.hpp>

namespace sgpp {
namespace base {

double OperationQuadratureModLinear::doQuadrature(DataVector& alpha) {
  double res = 0;

  double tmp;
  int nr_outermost;  // nr. of outermost ("folded up") basis functions
  int cur_ind, cur_lev;
  GridStorage::point_type index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; ++iter) {
    tmp = pow(2.0, -static_cast<double>(iter->first->getLevelSum())) * alpha.get(iter->second);
    nr_outermost = 0;

    for (size_t d = 0; d < iter->first->getDimension(); ++d) {
      cur_ind = iter->first->getIndex(d);
      cur_lev = iter->first->getLevel(d);

      // first and last index have outermost ("folded up") basis functions => double the area
      if ((cur_ind == 1) || (pow(2.0, cur_lev) - 1 == cur_ind)) ++nr_outermost;
    }

    tmp *= (pow(2.0, nr_outermost));

    res += tmp;
  }

  /*
    double tmpres = 0;

    for (size_t i = 0; i < alpha.getSize(); i++) {
      GridPoint& gp = storage.getPoint(i);
      tmpres = 1.;

      for (size_t d = 0; d < storage.getDimension(); d++) {
        tmpres *= base.getIntegral(gp.getLevel(d), gp.getIndex(d));
      }
      res += alpha[i] * tmpres;
    }
  */

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  for (size_t d = 0; d < storage.getDimension(); d++) {
    res *= storage.getBoundingBox()->getIntervalWidth(d);
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
