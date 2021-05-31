// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Class for representing a collection of full grids together with one scalar coefficient
 * per full grid.
 */
class CombinationGrid {
 public:
  /**
   * Default constructor, leaves the list of full grids empty.
   */
  CombinationGrid();

  /**
   * Constructor.
   *
   * @param fullGrids     vector of full grids
   * @param coefficients  vector of coefficients, same size as \c fullGrids
   */
  CombinationGrid(const std::vector<FullGrid>& fullGrids, const base::DataVector& coefficients);

  /**
   * Constructor for a CombinationGrid consisting of only one FullGrid.
   *
   * @param fullGrid      the full grid
   */
  explicit CombinationGrid(const FullGrid& fullGrid);

  /**
   * Factory method to create a CombinationGrid corresponding to the combination technique for
   * a regular sparse grid.
   *
   * @param dim           dimensionality
   * @param n             sparse grid level
   * @param basis         basis of the sparse grid (will be the same for all full grids;
   *                      this can be changed)
   * @param hasBoundary   whether the sparse grid has points on the boundary
   * @return CombinationGrid corresponding to the combination technique for the regular sparse grid
   */
  static CombinationGrid fromRegularSparse(size_t dim, level_t n, const HeterogeneousBasis& basis,
                                           bool hasBoundary = true);

  /**
   * Factory method to create a CombinationGrid corresponding to the truncated combination technique
   * for a regular sparse grid, given a truncation level
   *
   * @param dim           dimensionality
   * @param truncationLevel minimum level for each full component grid
   * @param levelSumDistance regular level in analogy to `fromRegularSparse`, the result is
   *                      a regular sparse grid of level `n` with each level shifted by
   *                      `truncationLevel`
   * @param basis         basis of the sparse grid (will be the same for all full grids;
   *                      this can be changed)
   * @param hasBoundary   whether the sparse grid has points on the boundary
   * @return CombinationGrid corresponding to the combination technique for the regular sparse grid
   */
  static CombinationGrid fromRegularSparseTruncated(size_t dim, LevelVector truncationLevel,
                                                    level_t levelSumDistance,
                                                    const HeterogeneousBasis& basis,
                                                    bool hasBoundary = true);

  /**
   * Factory method to create a CombinationGrid corresponding to the combination technique for
   * a (potentially dimensionally adaptive) sparse grid specified via its
   * nodal subspaces/full grids.
   *
   * Note: You have to specify the levels of *all* nodal subspaces, i.e., the specified list of
   * subspaces has to be downward closed.
   *
   * @param subspaceLevels  vector of subspace levels
   * @param basis           basis of the sparse grid (will be the same for all full grids;
   *                        this can be changed)
   * @param hasBoundary     whether the sparse grid has points on the boundary
   * @return CombinationGrid corresponding to the combination technique for the sparse grid
   */
  static CombinationGrid fromSubspaces(const std::vector<LevelVector>& subspaceLevels,
                                       const HeterogeneousBasis& basis, bool hasBoundary = true);

  /**
   * Combine the grid points of all full grids of this combination grid and store the grid points
   * in an existing GridStorage.
   *
   * @param[out] gridStorage  GridStorage in which to save the grid points (will be cleared)
   */
  void combinePoints(base::GridStorage& gridStorage) const;

  /**
   * Combine scalars (one scalar per full grid) using a weighted sum
   * (weighted by the coefficients of the combination grid).
   *
   * @param values  data vector, same size as the number of full grids
   * @return value resulting from the combination
   */
  double combineValues(const base::DataVector& values) const;

  /**
   * Combine equally-sized vectors (one per full grid) using a weighted sum
   * (weighted by the coefficients of the combination grid).
   *
   * @param[in] values    data matrix, same number of columns as the number of full grids
   *                      (every column corresponds to one full grid, the order of rows is given by
   *                      IndexVectorRange)
   * @param[out] result   vector resulting from the combination
   */
  void combineValues(const base::DataMatrix& values, base::DataVector& result) const;

  /**
   * Combine scalars associated to every full grid point using a weighted sum
   * (weighted by the coefficients of the combination grid).
   *
   * Hierarchization is one example: In every full grid, there is a scalar (hierarchical surplus)
   * associated with every grid point. This method combines these grid values to have one value
   * for each point in the combined grid (usually a sparse grid).
   *
   * @param[in] gridStorage   GridStorage containing the combined grid
   * @param[in] values        vector of DataVector, each DataVector corresponds to a full grid
   *                          and has the same size as the number of grid points of the full grid
   *                          (the order of DataVector entries is given by IndexVectorRange)
   * @param[out] result       vector resulting from the combination, same order as \c gridStorage
   */
  void combineSparseGridValues(const base::GridStorage& gridStorage,
                               const std::vector<base::DataVector>& values,
                               base::DataVector& result) const;

  /**
   * Combine vectors associated to every full grid point using a weighted sum
   * (weighted by the coefficients of the combination grid), i.e., vector version of the other
   * \c combineSparseGridValues method.
   *
   * @param[in] gridStorage   GridStorage containing the combined grid
   * @param[in] values        vector of DataMatrix, each DataMatrix corresponds to a full grid
   *                          and has the same number of columns as the number of grid points of
   *                          the full grid (every column corresponds to one full grid point,
   *                          the order of columns is given by IndexVectorRange)
   * @param[out] result       matrix resulting from the combination, columns have the same order
   *                          as \c gridStorage (every column corresponds to one grid point of the
   *                          combined grid)
   */
  void combineSparseGridValues(const base::GridStorage& gridStorage,
                               const std::vector<base::DataMatrix>& values,
                               base::DataMatrix& result) const;

  /**
   * Distribute values given on the combined grid to a specific full grid
   * (which should be contained in this combination grid).
   *
   * @param[in] gridStorage   GridStorage containing the combined grid
   * @param[in] values        vector of values on the combined grid, same size as \c gridStorage
   * @param[in] fullGrid      full grid, should be contained in this combination grid
   * @param[out] result       vector of values on the full grid, same size as the number of grid
   *                          points of the full grid (the order is given by IndexVectorRange)
   */
  void distributeValuesToFullGrid(const base::GridStorage& gridStorage,
                                  const base::DataVector& values, const FullGrid& fullGrid,
                                  base::DataVector& result) const;

  /**
   * Distribute values given on the combined grid to the full grids contained in this combination
   * grid.
   *
   * @param[in] gridStorage   GridStorage containing the combined grid
   * @param[in] values        vector of values on the combined grid, same size as \c gridStorage
   * @param[out] result       vector of vectors with values on the full grids, every vector
   *                          corresponds to one full grid of the combination grid, every vector
   *                          has the same size as the number of grid points of the respective
   *                          full grid (the order of DataVector entries is given by
   *                          IndexVectorRange)
   */
  void distributeValuesToFullGrids(const base::GridStorage& gridStorage,
                                   const base::DataVector& values,
                                   std::vector<base::DataVector>& result) const;

  /**
   * @return dimensionality
   */
  size_t getDimension() const;

  /**
   * @return vector of full grids
   */
  const std::vector<FullGrid>& getFullGrids() const;

  /**
   * @return vector of coefficients, same size as \c fullGrids
   */
  const base::DataVector& getCoefficients() const;

  /**
   * @param fullGrids vector of full grids
   * @param coefficients vector of coefficients, same size as \c fullGrids
   */
  void setFullGridsAndCoefficients(const std::vector<FullGrid>& fullGrids,
                                   const base::DataVector& coefficients);

 protected:
  /// vector of full grids
  std::vector<FullGrid> fullGrids;
  /// vector of coefficients, same size as \c fullGrids
  base::DataVector coefficients;
};

base::DataVector getStandardCoefficientsFromLevelSet(const std::vector<LevelVector>& levelSet);

}  // namespace combigrid
}  // namespace sgpp
