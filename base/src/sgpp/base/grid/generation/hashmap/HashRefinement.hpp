// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * Free refinement class for sparse grids
 */
class HashRefinement: public AbstractRefinement {
 public:
  /**
   * Refines a grid according to a RefinementFunctor provided.
   * Refines up to RefinementFunctor::getRefinementsNum() grid points if
   * possible, and if their refinement value is larger than RefinementFunctor::start()
   * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
   * If addedPoints is supplied and not a zero pointer, then newly created grid points are
   * appended to this vector.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a RefinementFunctor specifying the refinement criteria
   * @param addedPoints pointer to vector to append newly created grid points to
   */
  void free_refine(GridStorage& storage,
                   RefinementFunctor& functor,
                   std::vector<size_t>* addedPoints = nullptr) override;


  /**
   * Computes and returns the number of grid points, which can be refined.
   * This is the number of grid points that have at least one child missing.
   *
   * @param storage hashmap that stores the grid points
   * @return The number of grid points that can be refined
   */
  size_t getNumberOfRefinablePoints(GridStorage& storage) override;

  /**
   * Refine one grid point along a single direction
   * @param storage hashmap that stores the grid points
   * @param point point to refine
   * @param d direction
   */
  void refineGridpoint1D(GridStorage& storage, GridPoint& point, size_t d) override;
  void refineGridpoint1D(GridStorage& storage, size_t seq, size_t d) override;

  ~HashRefinement() override {}


 protected:
  /**
   * This method refines a grid point by generating the children in every dimension
   * of the grid and all their missing ancestors by calling create_gridpoint().
   *
   * @param storage hashmap that stores the gridpoints
   * @param refine_index The index in the hashmap of the point that should be refined
   */
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;

  /**
   * This method creates a new point on the grid. It checks if some parents or
   * children are needed in other dimensions.
   *
   * @param storage hashmap that stores the gridpoints
   * @param point The point that should be inserted
   */
  void createGridpoint(GridStorage& storage, GridPoint& point) override;

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

#endif /* HASHREFINEMENT_HPP */
