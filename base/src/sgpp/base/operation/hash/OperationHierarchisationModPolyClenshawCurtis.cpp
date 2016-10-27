// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationModPolyClenshawCurtis.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationModPolyClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationModPolyClenshawCurtis.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationModPolyClenshawCurtis::doHierarchisation(DataVector& node_values) {
  HierarchisationModPolyClenshawCurtis func(storage, &base);
  sweep<HierarchisationModPolyClenshawCurtis> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationModPolyClenshawCurtis::doDehierarchisation(DataVector& alpha) {
  DehierarchisationModPolyClenshawCurtis func(storage, &base);
  sweep<DehierarchisationModPolyClenshawCurtis> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    DataVector source(alpha);
    s.sweep1D(source, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
