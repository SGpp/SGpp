// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>

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
 * Takes the levels from levelStructure and stores its point in the given HashGridStorage (points
 * with boundary)
 *
 * @param levelStructure level structure of anisotropic full grids of combination technique
 * @param storage hash map of sparse grid points
 */
void convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(
    std::shared_ptr<TreeStorage<uint8_t>> levelStructure, base::HashGridStorage& storage);

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

/**
 * converts a given level structure from the combigrid module into a matrix containing all grid
 * points from these levels
 * @param levelStructure the level structure
 * @param numDimensions number of dimensions
 * @return matrix containing the grid points
 */
sgpp::base::DataMatrix convertLevelStructureToGridPoints(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    size_t numDimensions);

/**
 * after a combigrid has been converted to a hierarchical sparse grid via
 * convertCombigridToHierarchicalSparseGrid with this operation the coefficients for the
 * hierarchical sparse grid interpolation can be calculated
 * @param grid the newly created hierarchical sparse grid
 * @param gridStorage the gridStorage of the newly created hierarchical sparse grid
 * @param combigridInterpolationOperation an interpolation operation on the comigrid used to obtain
 * the right hand side for the interpolation SLE
 * @param levelStructure the level structure of the combigrid
 * @return the coefficients for the hierarchical sparse grid interpolant
 */
sgpp::base::DataVector calculateInterpolationCoefficientsForConvertedExpUniformBoundaryCombigird(
    std::shared_ptr<sgpp::base::Grid>& grid, sgpp::base::GridStorage& gridStorage,
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation>& combigridInterpolationOperation,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure);

} /* namespace combigrid */
} /* namespace sgpp */
