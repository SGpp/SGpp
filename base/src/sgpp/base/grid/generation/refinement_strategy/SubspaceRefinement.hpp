// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SUBSPACEREFINEMENT_HPP_
#define SUBSPACEREFINEMENT_HPP_

#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <iostream>


namespace sgpp {
namespace base {


/*
 * SubspaceRefinement inserts complete hierarchical subspaces to the sparse
 * grid instead of individual grid points.
 */
class SubspaceRefinement: public RefinementDecorator {
 public:
  /**
   * Constructor
   *
   * @param refinement decorated refinement object
   */
  explicit SubspaceRefinement(AbstractRefinement* refinement):
    RefinementDecorator(refinement) {
  }


  /**
  * Refines a grid according to a RefinementFunctor provided.
  * Refines up to RefinementFunctor::getRefinementsNum() grid points if
  * possible, and if their refinement value is larger than RefinementFunctor::start()
  * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
  *
  * @param storage hashmap that stores the grid points
  * @param functor a RefinementFunctor specifying the refinement criteria
  * @param addedPoints pointer to vector to append newly created grid points to
  */
  void free_refine(GridStorage& storage,
                   RefinementFunctor& functor,
                   std::vector<size_t>* addedPoints = nullptr) override;


  /**
   * Destructor
   */
  ~SubspaceRefinement() override {}

 protected:
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
  * @param refinement_num number of elements to refine
  * @param collection container where element pairs for refinement need to be stored
  */
  virtual void addElementToCollection(
    const GridStorage::grid_map_iterator& iter,
    AbstractRefinement::refinement_list_type current_value_list,
    size_t refinement_num,
    AbstractRefinement::refinement_container_type& collection);
};

}  // namespace base
}  // namespace sgpp
#endif /* SUBSPACEREFINEMENT_HPP_ */
