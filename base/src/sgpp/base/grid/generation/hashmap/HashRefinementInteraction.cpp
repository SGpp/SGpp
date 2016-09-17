// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInteraction.hpp>

#include <algorithm>
#include <memory>
#include <vector>

namespace sgpp {
namespace base {

HashRefinementInteraction::HashRefinementInteraction(
    std::unordered_set<std::vector<bool>> interactions)
    : HashRefinement(), interactions(interactions) {}

void HashRefinementInteraction::createGridpoint(GridStorage& storage, index_type& index) {
  index_t source_index;
  level_t source_level;

  // Get the currently used dimensions of the grid point.
  auto coordsIndex = std::vector<bool>(storage.getDimension());
  for (size_t i = 0; i < storage.getDimension(); ++i) {
    coordsIndex[i] = index.getCoord(i) != 0.5;
  }

  // Check first whether index obeys our limitations.
  if (interactions.find(coordsIndex) == interactions.end()) {
    return;
  }

  // We now check whether the index with additional dimension d is legal.
  for (size_t d = 0; d < storage.getDimension(); d++) {
    auto curCoords = coordsIndex;
    curCoords[d] = true;
    if (interactions.find(curCoords) != interactions.end()) {
      createGridpoint1D(index, d, storage, source_index, source_level);
    }
  }

  storage.insert(index);
}

}  // namespace base
}  // namespace sgpp
