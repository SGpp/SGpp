// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDotProductModLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationDotProductModLinear::eval(base::DataVector& x1, base::DataVector& x2) {
  base::LinearModifiedBasis<unsigned int, unsigned int> base;
  base::GridStorage::index_type::level_type work_level = 1;
  base::GridStorage::index_type::index_type work_index;
  base::GridStorage::index_type::level_type temp;
  double result = 0;

  // GridStorage::grid_iterator working;
  // for (GridStorage::grid_iterator working = storage->begin(); working != storage->end();
  // working++){
  for (size_t i = 0; i < storage->getSize(); i++) {
    base::GridStorage::index_type& working = storage->getGridIndex(i);
    double value1 = 1.0;
    double value2 = 1.0;

    for (size_t d = 0; d < storage->getDimension(); d++) {
      working.get(d, temp, work_index);

      value1 *= base.eval(work_level, work_index, x1[d]);
      value2 *= base.eval(work_level, work_index, x2[d]);

      //}
    }

    result += value1 * value2;
  }

  return result;
}
}  // namespace datadriven
}  // namespace sgpp
