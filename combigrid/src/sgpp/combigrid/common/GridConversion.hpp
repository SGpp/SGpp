// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Creates a level structure from all levels that occur in the HashGridStorage.
 */
std::shared_ptr<TreeStorage<uint8_t>> allStorageLevels(base::HashGridStorage* storage);

/**
 * Creates a level structure from all levels of which all points occur in the HashGridStorage.
 * This assumes that at a full level (l_1, ..., l_n), there are 2^(l_1 + ... + l_n) points in the
 * HashGridStorage.
 */
std::shared_ptr<TreeStorage<uint8_t>> fullStorageLevels(base::HashGridStorage* storage);

/**
 * Takes the levels from levelStructure and stores its point in the given HashGridStorage (points
 * without boundary).
 */
void toHashGridStorage(std::shared_ptr<TreeStorage<uint8_t>> levelStructure,
                       base::HashGridStorage& storage);

} /* namespace combigrid */
} /* namespace sgpp */
