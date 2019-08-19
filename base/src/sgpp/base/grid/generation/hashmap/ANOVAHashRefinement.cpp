// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/ANOVAHashRefinement.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace base {

void ANOVAHashRefinement::refineGridpoint(GridStorage& storage, size_t refine_index) {
  GridPoint point(storage[refine_index]);
  // Sets leaf property of index, which is refined to false
  storage[refine_index].setLeaf(false);

  for (size_t d = 0; d < storage.getDimension(); d++) {
    // For ANOVA refinement create children only in the dimensions with level
    // greater than 1 (non-constant basis functions)
    if (point.getLevel(d) > 1) {
      this->refineGridpoint1D(storage, point, d);
    }
  }
}
}  // namespace base
}  // namespace sgpp
