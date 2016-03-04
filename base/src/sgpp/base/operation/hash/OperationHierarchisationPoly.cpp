// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationPoly.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPoly.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationPoly.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationPoly::doHierarchisation(DataVector& node_values) {
  HierarchisationPoly func(storage, &base);
  sweep<HierarchisationPoly> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationPoly::doDehierarchisation(DataVector& alpha) {
  DehierarchisationPoly func(storage, &base);
  sweep<DehierarchisationPoly> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    DataVector source(alpha);
    s.sweep1D(source, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
