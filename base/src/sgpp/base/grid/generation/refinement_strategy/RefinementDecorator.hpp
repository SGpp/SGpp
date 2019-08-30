// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REFINEMENTSTRATEGY_HPP_
#define REFINEMENTSTRATEGY_HPP_

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * RefinementDecorator enhances the behavior of underlying Refinement objects
 * using <a href="http://en.wikipedia.org/wiki/Decorator_pattern"> Decorator design
 * pattern </a>. Although not abstract, this class is thought
 * to be a base class as it simply delegates the function calls to the decorated
 * object. Subclasses will implement more sophisticated behavior.
 */
class RefinementDecorator: public AbstractRefinement {
 public:
  /**
   * Constructor
   *
   * @param refinement object implementing the core functionality (e.g.
   * refinement with or without boundaries).
   */
  explicit RefinementDecorator(AbstractRefinement* refinement) {
    decorated_refinement_ = refinement;
  }

  virtual ~RefinementDecorator() {
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
  virtual void free_refine(GridStorage& storage, RefinementFunctor& functor,
                           std::vector<size_t>* addedPoints = nullptr);


  /**
   * Computes and returns the number of grid points, which can be refined.
   * This is the number of grid points that have at least one child missing.
   *
   * @param storage hashmap that stores the grid points
   * @return The number of grid points that can be refined
   */
  virtual size_t getNumberOfRefinablePoints(GridStorage& storage);

  /**
   * Refine one grid point along a single direction
   * @param storage hashmap that stores the grid points
   * @param point point to refine
   * @param d direction
   */
  virtual void refineGridpoint1D(GridStorage& storage, GridPoint& point,
                                 size_t d);

  bool checkAdmissibility(GridStorage& storage, GridPoint& subspace);

 protected:
  /**
   * Returns the pointer to decorated Refinement object
   */
  AbstractRefinement* get_decorated_refinement() {
    return decorated_refinement_;
  }

  /**
   * Sets the pointer of the decorated Refinement object to the given object
   *
   * @param refinement object the pointer should be set to
   */
  void set_decorated_refiment(AbstractRefinement* refinement) {
    decorated_refinement_ = refinement;
  }

  /**
   * This method refines a grid point by generating the children in every dimension
   * of the grid and all their missing ancestors by calling create_gridpoint().
   *
   * @param storage hashmap that stores the gridpoints
   * @param refine_index The index in the hashmap of the point that should be refined
   */
  virtual void refineGridpoint(GridStorage& storage, size_t refine_index);

  /**
   * This method creates a new point on the grid. It checks if some parents or
   * children are needed in other dimensions.
   *
   * @param storage hashmap that stores the gridpoints
   * @param point The point that should be inserted
   */
  virtual void createGridpoint(GridStorage& storage, GridPoint& point);

  /**
  * Examines the grid points and stores the indices those that can be refined
  * and have maximal indicator values.
  *
  * @param storage hashmap that stores the grid points
  * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
  * @param collection container that contains elements to refine (empty initially)
  */
  virtual void collectRefinablePoints(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection);


  /**
   * Refines the collection of points.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
   * @param collection container that contains elements to refine (empty initially)
   */
  virtual void refineGridpointsCollection(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection);



  /**
   * Generates a list with indicator elements
   *
   * @param storage grid storage
   * @param iter iterator
   * @param functor refinement functor
   * @return list with indicator elements
   */
  virtual AbstractRefinement::refinement_list_type getIndicator(
    GridStorage& storage,
    const GridStorage::grid_map_iterator& iter,
    const RefinementFunctor& functor) const;

 private:
  AbstractRefinement* decorated_refinement_;
};

}  // namespace base
}  // namespace sgpp

#endif /* REFINEMENTSTRATEGY_HPP_ */
