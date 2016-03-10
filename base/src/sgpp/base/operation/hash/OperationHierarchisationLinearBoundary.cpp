// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationHierarchisationLinearBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLinearBoundary func(storage);
  sweep<HierarchisationLinearBoundary> s(func, storage);

  // N D case
  if (this->storage.getDimension() > 1) {
    for (size_t i = 0; i < this->storage.getDimension(); i++) {
      s.sweep1D_Boundary(node_values, node_values, i);
    }
  } else {  // 1 D case
    s.sweep1D(node_values, node_values, 0);
  }
}

void OperationHierarchisationLinearBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationLinearBoundary func(storage);
  sweep<DehierarchisationLinearBoundary> s(func, storage);

  // N D case
  if (this->storage.getDimension() > 1) {
    for (size_t i = 0; i < this->storage.getDimension(); i++) {
      s.sweep1D_Boundary(alpha, alpha, i);
    }
  } else {  // 1 D case
    s.sweep1D(alpha, alpha, 0);
  }
}

}  // namespace base
}  // namespace sgpp
