// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void OperationHierarchisationLinearBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLinearBoundary func(this->storage);
  sweep<HierarchisationLinearBoundary> s(func, this->storage);

  // N D case
  if (this->storage->dim() > 1) {
    for (size_t i = 0; i < this->storage->dim(); i++) {
      s.sweep1D_Boundary(node_values, node_values, i);
    }
  }
  // 1 D case
  else {
    s.sweep1D(node_values, node_values, 0);
  }
}

void OperationHierarchisationLinearBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationLinearBoundary func(this->storage);
  sweep<DehierarchisationLinearBoundary> s(func, this->storage);

  // N D case
  if (this->storage->dim() > 1) {
    for (size_t i = 0; i < this->storage->dim(); i++) {
      s.sweep1D_Boundary(alpha, alpha, i);
    }
  }
  // 1 D case
  else {
    s.sweep1D(alpha, alpha, 0);
  }
}

}
}