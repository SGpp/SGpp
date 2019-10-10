// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ONLINEPREDICTIVEREFINEMENTDIMENSIONOLD_HPP_
#define ONLINEPREDICTIVEREFINEMENTDIMENSIONOLD_HPP_

#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <vector>
#include <utility>


namespace sgpp {
namespace base {

/**
 * Container type for predictive refinement collection
 */
class PredictiveRefinement_refinement_key : public
  AbstractRefinement_refinement_key {
 public:
  /**
   * Constructor
   *
   * @param point grid point
   * @param seq sequence number in the hash grid storage
   * @param dim dimensionality
   */
  PredictiveRefinement_refinement_key(const GridPoint& point, size_t seq, size_t dim):
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
  virtual ~PredictiveRefinement_refinement_key() {}

 private:
  size_t dim;
};


/*
 * PredictiveRefinement performs local adaptive refinement of a sparse grid using the PredictiveRefinementIndicator.
 * This way, instead of surpluses that are most often used for refinement, refinement decisions are based upon an estimation
 * to the contribution of the MSE, which is especially helpful for regression.
 *
 */
class PredictiveRefinement: public virtual RefinementDecorator {
  friend class LearnerOnlineSGD;
 public:
  typedef PredictiveRefinement_refinement_key refinement_key_type;
  using RefinementDecorator::free_refine;

  explicit PredictiveRefinement(AbstractRefinement* refinement):
    RefinementDecorator(refinement),
    iThreshold_(0.0), alpha_(0) {}


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
                   PredictiveRefinementIndicator& functor,
                   std::vector<size_t>* addedPoints = nullptr);





  /**
   * Setter for the alpha vector
   * @param alpha
   */
  void setAlpha(DataVector& alpha) {
    alpha_ = alpha;
  }

 protected:
  using RefinementDecorator::refineGridpointsCollection;

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
    AbstractRefinement::refinement_container_type&  collection) override;


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
  DataVector alpha_;
};

}  // namespace base
}  // namespace sgpp
#endif /* ONLINEPREDICTIVEREFINEMENTDIMENSIONOLD_HPP_ */
