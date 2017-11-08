// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

namespace sgpp {
namespace combigrid {

enum class GridConversionTypes { ALLSUBSPACES, COMPLETESUBSPACES };

/**
 * Creates a level structure for the levels that occur in the HashGridStorage.
 * Depending on the conversion type, it considers all available subspaces or
 * just the complete ones.
 *
 * @param storage
 * @param conversionType
 * @return level structure of anisotropic full grids of combination technique
 */
std::shared_ptr<TreeStorage<uint8_t>> convertHierarchicalSparseGridToCombigrid(
    base::HashGridStorage& storage, GridConversionTypes conversionType);

/**
 * Takes the levels from levelStructure and stores its point in the given HashGridStorage (points
 * without boundary)
 *
 * @param levelStructure level structure of anisotropic full grids of combination technique
 * @param storage hash map of sparse grid points
 */
void convertCombigridToHierarchicalSparseGrid(std::shared_ptr<TreeStorage<uint8_t>> levelStructure,
                                              base::HashGridStorage& storage);

/**
 * Creates a level structure from all levels that occur in the HashGridStorage.
 *
 * @param storage hash map of sparse grid points
 * @return level structure of anisotropic full grids for CombigridOperations
 */
std::shared_ptr<TreeStorage<uint8_t>> allStorageLevels(base::HashGridStorage& storage);

/**

 */
/**
 * Creates a level structure from all levels of which all points occur in the HashGridStorage.
 * This assumes that at a full level (l_1, ..., l_n), there are 2^(l_1 + ... + l_n) points in the
 * HashGridStorage.
 *
 * @param storage hash map of sparse grid points
 * @return level structure of anisotropic full grids for CombigridOperations
 */
std::shared_ptr<TreeStorage<uint8_t>> completeStorageLevels(base::HashGridStorage& storage);

} /* namespace combigrid */
} /* namespace sgpp */
