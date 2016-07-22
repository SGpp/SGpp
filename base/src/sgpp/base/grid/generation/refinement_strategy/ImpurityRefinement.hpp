// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef IMPURITYREFINEMENT_HPP_
#define IMPURITYREFINEMENT_HPP_

#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

#include <vector>
#include <utility>


namespace sgpp {
namespace base {

/**
 * Container type for impurity refinement collection
 */
class ImpurityRefinement_refinement_key : public
  AbstractRefinement_refinement_key {
 public:
  /**
   * Constructor
   *
   * @param point grid point
   * @param seq sequence number in the hash grid storage
   * @param dim dimensionality
   */
  ImpurityRefinement_refinement_key(const GridPoint& point, size_t seq, size_t dim):
    AbstractRefinement_refinement_key(point, seq), dim(dim) {}

  /**
   * Returns dimensionality
   *
   * @return dimensionality
   */
  size_t getDim() {
    return this->dim;
  }

  /**
   * Destructor
   */
  virtual ~ImpurityRefinement_refinement_key() {}

 private:
  size_t dim;
};


/*
 * 
 * 
 * 
 *
 */
class ImpurityRefinement: public virtual RefinementDecorator {
 public:
  typedef ImpurityRefinement_refinement_key refinement_key_type;
  using RefinementDecorator::free_refine;

  explicit ImpurityRefinement(AbstractRefinement* refinement):
    RefinementDecorator(refinement),
    iThreshold_(0.0) {}


  /**
   * Refines a grid according to a RefinementFunctor provided.
   * Refines up to RefinementFunctor::getRefinementsNum() grid points if
   * possible, and if their refinement value is larger than RefinementFunctor::start()
   * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
   *
   * @param storage hashmap that stores the grid points
   * @param functor a RefinementFunctor specifying the refinement criteria
   */
  void free_refine(GridStorage& storage,
                   ImpurityRefinementIndicator& functor);


 protected:
  using RefinementDecorator::refineGridpointsCollection;

  /**
  * Examines the grid points and stores the indices those that can be refined
  * and have maximal indicator values.
  *
  * @param storage hashmap that stores the grid points
  * @param functor an impurity indicator specifying the refinement criteria
  * @param collection container that contains elements to refine (empty initially)
  */
  void collectRefinablePoints(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type&  collection) override;


  /**
   * Extends the grid adding elements defined in collection
   *
   * @param storage hashmap that stores the grid points
   * @param functor an impurity indicator specifying the refinement criteria
   * @param collection container that contains elements to refine (empty initially)
   */
  void refineGridpointsCollection(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) override;

  void createGridpointMod(GridStorage& storage, GridPoint& point, 
    ImpurityRefinementIndicator& indicator);

  void createGridpoint1DMod(GridPoint& point, size_t d, GridStorage& storage,
    index_t& source_index, level_t& source_level,
    ImpurityRefinementIndicator& indicator);  

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

 private:
  double iThreshold_;
};

}  // namespace base
}  // namespace sgpp
#endif /* IMPURITYREFINEMENT_HPP_ */
