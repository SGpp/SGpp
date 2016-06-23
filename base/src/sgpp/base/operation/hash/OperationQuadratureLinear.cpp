// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureLinear.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

double OperationQuadratureLinear::doQuadrature(DataVector& alpha) {
  double res = 0;
  GridStorage::point_type index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    //    index = *(iter->first);
    //    std::cout << iter->second << " " << iter->first->getLevelSum() <<
    //    " " << pow(2.0, -static_cast<double>(iter->first->getLevelSum())) <<
    //    std::endl;
    res += pow(2.0, -static_cast<double>(iter->first->getLevelSum())) * alpha.get(
             iter->second);
  }

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  for (size_t d = 0; d < storage.getDimension(); d++) {
    res *= storage.getBoundingBox()->getIntervalWidth(d);
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
