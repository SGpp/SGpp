// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinearClenshawCurtis.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearClenshawCurtis.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationLinearClenshawCurtis::doHierarchisation(DataVector& node_values) {
  HierarchisationLinearClenshawCurtis func(storage);
  sweep<HierarchisationLinearClenshawCurtis> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationLinearClenshawCurtis::doDehierarchisation(DataVector& alpha) {
  DehierarchisationLinearClenshawCurtis func(storage);
  sweep<DehierarchisationLinearClenshawCurtis> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    DataVector source(alpha);
    s.sweep1D(source, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
