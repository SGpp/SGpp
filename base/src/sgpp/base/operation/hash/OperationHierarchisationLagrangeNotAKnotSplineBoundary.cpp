// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLagrangeNotAKnotSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLagrangeNotAKnotSplineBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLagrangeNotAKnotSplineBoundary func(storage);
  sweep<HierarchisationLagrangeNotAKnotSplineBoundary> s(func, storage);

  // N D case
  if (this->storage.getDimension() > 1) {
    for (size_t i = 0; i < this->storage.getDimension(); i++) {
      s.sweep1D_Boundary(node_values, node_values, i);
    }
  } else {  // 1 D case
    s.sweep1D(node_values, node_values, 0);
  }
}

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  throw sgpp::base::not_implemented_exception("Dehierarchisation is not yet implemented!");
}

}  // namespace base
}  // namespace sgpp
