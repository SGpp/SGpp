// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENTBOUNDARIES_HPP
#define HASHREFINEMENTBOUNDARIES_HPP

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * Standard free refinement class for sparse grids with boundaries
 */
class HashRefinementBoundaries: public AbstractRefinement {
 public:
  /**
   * Performs the refinement on grid
   *
   * @param storage hashmap that stores the grid points
   * @param functor a function used to determine if refinement is needed
   * @param addedPoints pointer to vector to append newly created grid points to
   */
  void free_refine(GridStorage& storage,
                   RefinementFunctor& functor,
                   std::vector<size_t>* addedPoints = nullptr) override;


  /**
   * Calculates the number of points, which can be refined
   *
   * @param storage hashmap that stores the grid points
   */
  size_t getNumberOfRefinablePoints(GridStorage& storage) override;

  void refineGridpoint1D(GridStorage& storage, GridPoint& point, size_t d) override;

 protected:
  /**
   * This method refines a grid point be generating the children in every dimension
   * of the grid.
   *
   * @param storage hashmap that stores the gridpoints
   * @param refine_index the index in the hashmap of the point that should be refined
   */
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;


  /**
   * Wrapper for the two functions create_gridpoint_general and
   * create_gridpoint_levelZeroConsistency which have both to be
   * executed if a gridpoint is refined
   *
   * @param storage hashmap that stores the gridpoinrs
   * @param point the point that should be inserted
   */
  void createGridpoint(GridStorage& storage, GridPoint& point) override;


  /**
   * This method creates a new point on the grid. It checks if some parents or
   * children are needed in other dimensions.
   *
   * @param storage hashmap that stores the gridpoinrs
   * @param point the point that should be inserted
   */
  void createGridpointGeneral(GridStorage& storage, GridPoint& point);


  /**
   * Assures that we have always a consistent grid with both functions
   * 0,0 and 0,1 on level zero
   *
   * @param storage hashmap that stores the gridpoinrs
   * @param point the point that should be inserted
   */
  void createGridpointLevelZeroConsistency(GridStorage& storage, GridPoint& point);

  /**
         * Creates children grid points along single direction
         *
         * @param point The point that should be refined
         * @param d direction
         * @param storage hashmap that stores the gridpoints
         * @param source_index index value in the dimension d
         * @param source_level level value in the dimension d
         */
  void createGridpoint1D(GridPoint& point,
                         size_t d, GridStorage& storage,
                         index_t& source_index, level_t& source_level) override;

  /**
  * Examines the grid points and stores the indices those that can be refined
  * and have maximal indicator values.
  *
  * @param storage hashmap that stores the grid points
  * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
  * @param collection container that contains elements to refine (empty initially)
  */
  void collectRefinablePoints(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) override;

  /**
  * Extends the grid adding elements defined in collection
  *
  * @param storage hashmap that stores the grid points
  * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
  * @param collection container that contains elements to refine (empty initially)
  */
  void refineGridpointsCollection(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) override;

  /**
  * Adds elements to the collection. This method is responsible for selection
  * the elements with most important indicators and to limit the size of collection
  * to refinements_num elements.
  *
  * @param iter storage iterator
  * @param current_value_list list with elements that contain keys and values that specify refinement
  * @param refinements_num number of elements to refine
  * @param collection container where element pairs for refinement need to be stored
  */
  virtual void addElementToCollection(
    const GridStorage::grid_map_iterator& iter,
    AbstractRefinement::refinement_list_type current_value_list,
    size_t refinements_num,
    AbstractRefinement::refinement_container_type& collection);

  /**
  * Generates a list with indicator elements
  *
  * @param storage grid storage
  * @param iter iterator
  * @param functor refinement functor
  * @return list with indicator elements
  */
  AbstractRefinement::refinement_list_type getIndicator(
    GridStorage& storage,
    const GridStorage::grid_map_iterator& iter,
    const RefinementFunctor& functor) const override;
};

}  // namespace base
}  // namespace sgpp

#endif /* HASHREFINEMENTBOUNDARIES_HPP */
