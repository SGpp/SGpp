// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHCOARSENING_HPP
#define HASHCOARSENING_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace base {

/**
 * Standard free coarsening class for sparse grids, only
 * inner grid points can be removed
 */
class HashCoarsening {
 public:
  /**
   * Performs coarsening on grid. It's possible to remove a certain number
   * of gridpoints in one coarsening step. This number is specified within the
   * declaration of the coarsening functor. Also the coarsening threshold is
   * specified in the coarsening functor. ONLY INNER GRID POINTS WILL
   * BE REMOVED!
   *
   * Here only the numFirstPoints are regarded for coarsening, later points
   * are skipped.
   *
   * Note that seq numbers in removedSeq are invalid to use with the resulting grid after
   * coarsening since the seq numbers get recomputed after removal. Use removedPoints unless
   * outdated seq numbers are explicitly required.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a function used to determine if refinement is needed
   * @param numFirstPoints number of grid points that are regarded to be coarsened
   * @param minIndexConsidered indices of coarsen point candidates must be higher than this
   * parameter to be allowed to get coarsened
   * @param removedPoints pointer to vector to append coarsened (removed) grid points to
   * @param removedSeq pointer to vector to append the seq numbers of coarsened grid points to
   */
  void free_coarsen_NFirstOnly(GridStorage& storage, CoarseningFunctor& functor,
                               size_t numFirstPoints, size_t minIndexConsidered = 0,
                               std::vector<HashGridPoint>* removedPoints = nullptr,
                               std::vector<size_t>* removedSeq = nullptr);

  /**
   * Performs coarsening on grid. It's possible to remove a certain number
   * of gridpoints in one coarsening step. This number is specified within the
   * declaration of the coarsening functor. Also the coarsening threshold is
   * specified in the coarsening functor. ONLY INNER GRID POINTS WILL
   * BE REMOVED!
   *
   * This function calls free_coarsen_NFirstOnly with numFirstPoints equal
   * to the grid's size.
   *
   * Note that seq numbers in removedSeq are invalid to use with the resulting grid after
   * coarsening since the seq numbers get recomputed after removal. Use removedPoints unless
   * outdated seq numbers are explicitly required.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a function used to determine if refinement is needed
   * @param removedPoints pointer to vector to append coarsened (removed) grid points to
   * @param removedSeq pointer to vector to append the seq numbers of coarsened grid points to.
   */
  void free_coarsen(GridStorage& storage, CoarseningFunctor& functor,
                    std::vector<HashGridPoint>* removedPoints = nullptr,
                    std::vector<size_t>* removedSeq = nullptr);

  /**
   * Calculates the number of points, which can be refined
   *
   * @param storage hashmap that stores the grid points
   */
  size_t getNumberOfRemovablePoints(GridStorage& storage);
};

}  // namespace base
}  // namespace sgpp

#endif /* HASHCOARSENING_HPP */
