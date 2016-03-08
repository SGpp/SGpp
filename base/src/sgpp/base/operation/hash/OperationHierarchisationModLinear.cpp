// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationModLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {
/**
 * Implements the hierarchisation on a sprase grid with mod linear base functions
 *
 * @param node_values the functions values in the node base
 *
 */
void OperationHierarchisationModLinear::doHierarchisation(
  DataVector& node_values) {
  HierarchisationModLinear func(storage);
  sweep<HierarchisationModLinear> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

/**
 * Implements the dehierarchisation on a sprase grid with mod linear base functions
 *
 * @param alpha the coefficients of the sparse grid's base functions
 *
 */
void OperationHierarchisationModLinear::doDehierarchisation(DataVector& alpha) {
  DehierarchisationModLinear func(storage);
  sweep<DehierarchisationModLinear> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(alpha, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
