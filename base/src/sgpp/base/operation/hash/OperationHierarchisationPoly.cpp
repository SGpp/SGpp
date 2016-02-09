// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "OperationHierarchisationPoly.hpp"

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPoly.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationPoly.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

void OperationHierarchisationPoly::doHierarchisation(DataVector& node_values) {
  HierarchisationPoly func(this->storage, &this->base);
  sweep<HierarchisationPoly> s(func, this->storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage->dim(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationPoly::doDehierarchisation(DataVector& alpha) {
  DehierarchisationPoly func(this->storage, &this->base);
  sweep<DehierarchisationPoly> s(func, this->storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage->dim(); i++) {
    DataVector source(alpha);
    s.sweep1D(source, alpha, i);
  }
}

}  // namespace base
}  // namespace SGPP
