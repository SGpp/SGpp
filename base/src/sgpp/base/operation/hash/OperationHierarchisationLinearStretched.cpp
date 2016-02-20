// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretched.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinearStretched.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinearStretched.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void OperationHierarchisationLinearStretched::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLinearStretched func(storage);
  sweep<HierarchisationLinearStretched> s(func, &storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.dim(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationLinearStretched::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationLinearStretched func(storage);
  sweep<DehierarchisationLinearStretched> s(func, &storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.dim(); i++) {
    s.sweep1D(alpha, alpha, i);
  }
}

}  // namespace base
}  // namespace SGPP
