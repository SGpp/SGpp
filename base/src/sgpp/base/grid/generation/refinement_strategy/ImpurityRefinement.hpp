// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef IMPURITYREFINEMENT_HPP_
#define IMPURITYREFINEMENT_HPP_

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * Container type for impurity refinement collection
 */
class ImpurityRefinement_refinement_key
    : public AbstractRefinement_refinement_key {
 public:
  /**
   * Constructor
   *
   * @param point grid point
   * @param seq sequence number in the hash grid storage
   * @param dim dimensionality
   */
  ImpurityRefinement_refinement_key(const GridPoint& point, size_t seq,
                                    size_t dim)
      : AbstractRefinement_refinement_key(point, seq), dim(dim) {}

  /**
   * Returns dimensionality
   *
   * @return dimensionality
   */
  size_t getDim() { return this->dim; }

  /**
   * Destructor
   */
  virtual ~ImpurityRefinement_refinement_key() {}

 private:
  size_t dim;
};

class ImpurityRefinement : public virtual RefinementDecorator {
 public:
  typedef ImpurityRefinement_refinement_key refinement_key_type;
  using RefinementDecorator::free_refine;

  explicit ImpurityRefinement(AbstractRefinement* refinement)
      : RefinementDecorator(refinement) {}

  /**
   * Refines a grid according to the impurity refinement indicator provided.
   *
   * @param storage Hashmap that stores the grid points
   * @param functor A RefinementFunctor specifying the refinement criteria
   */
  void free_refine(GridStorage& storage, ImpurityRefinementIndicator& functor);

 protected:
  using RefinementDecorator::refineGridpointsCollection;

  /**
  * Examines the grid points and stores the indices of those that can be refined
  * and have maximal indicator values.
  *
  * @param storage Hashmap that stores the grid points
  * @param functor An impurity indicator specifying the refinement criteria
  * @param collection Container that contains elements to refine (empty
  * initially)
  */
  void collectRefinablePoints(
      GridStorage& storage, RefinementFunctor& functor,
      AbstractRefinement::refinement_container_type& collection) override;

  /**
   * Extends the grid adding points defined in the collection
   *
   * @param storage Hashmap that stores the grid points
   * @param functor An impurity indicator specifying the refinement criteria
   * @param collection Container that contains elements to refine (empty
   * initially)
   */
  void refineGridpointsCollection(
      GridStorage& storage, RefinementFunctor& functor,
      AbstractRefinement::refinement_container_type& collection) override;

  /**
  * Generates a list with indicator elements
  *
  * @param storage Grid storage
  * @param iter Iterator
  * @param functor Refinement functor
  * @return List with indicator elements
  */
  AbstractRefinement::refinement_list_type getIndicator(
      GridStorage& storage, const GridStorage::grid_map_iterator& iter,
      const RefinementFunctor& functor) const override;

  /**
  * Adds elements to the collection. This method is responsible for selection of
  * the elements with largest indicator values and to limit the size of
  * collection
  * to refinementsNum elements.
  *
  * @param iter Storage iterator
  * @param current_value_list List with elements that contain keys and values
  * that specify refinement
  * @param refinementsNum Number of grid points to refine
  * @param collection Container where element pairs for refinement need to be
  * stored
  */
  virtual void addElementToCollection(
      const GridStorage::grid_map_iterator& iter,
      AbstractRefinement::refinement_list_type current_value_list,
      size_t refinementsNum,
      AbstractRefinement::refinement_container_type& collection);
};

}  // namespace base
}  // namespace sgpp

#endif /* IMPURITYREFINEMENT_HPP_ */
