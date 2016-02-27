// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationSecondMomentLinear.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

double OperationSecondMomentLinear::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmpres = 1;
  GridStorage::index_type index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    //    index = *(iter->first);
    //    std::cout << iter->second << " " << iter->first->getLevelSum() <<
    //    " " << pow(2.0, -static_cast<double>(iter->first->getLevelSum())) <<
    //    std::endl;
    tmpres = 1.;

    for (size_t dim = 0; dim < storage.getDimension(); dim++)
      tmpres *= pow(8.0, -static_cast<double>(iter->first->getLevel(dim))) *
                (iter->first->getIndex(dim) * iter->first->getIndex(dim) + 1. / 6.);

    res += alpha.get(iter->second) * tmpres;
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
