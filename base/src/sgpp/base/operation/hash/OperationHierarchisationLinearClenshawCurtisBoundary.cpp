// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinearClenshawCurtisBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearClenshawCurtisBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearClenshawCurtisBoundary.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationLinearClenshawCurtisBoundary::doHierarchisation(
    DataVector& node_values) {
  HierarchisationLinearClenshawCurtisBoundary func(storage);
  sweep<HierarchisationLinearClenshawCurtisBoundary> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D_Boundary(node_values, node_values, i);
  }
}

void OperationHierarchisationLinearClenshawCurtisBoundary::doDehierarchisation(DataVector& alpha) {
  DehierarchisationLinearClenshawCurtisBoundary func(storage);
  sweep<DehierarchisationLinearClenshawCurtisBoundary> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    DataVector source(alpha);
    s.sweep1D_Boundary(source, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
