// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ABSTRACTREFINEMENT_HPP
#define ABSTRACTREFINEMENT_HPP

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

#include <forward_list>
#include <iosfwd>
#include <vector>
#include <unordered_map>
#include <utility>
#include <memory>


namespace sgpp {
namespace base {

/***
 * Utility class for keys of the refinement containers.
 */
class AbstractRefinement_refinement_key {
 public:
  /***
   * Constructor
   *
   * @param point The grid point in the HashGridStorage
   * @param seq The sequential number of the grid points in the HashGridStorage
   */
  AbstractRefinement_refinement_key(const GridPoint& point, size_t seq) :
    point(point), seq(seq), level_vector() {
  }


  /***
   * Destructor
   */
  virtual ~AbstractRefinement_refinement_key() {}


  /***
   * Gets the level vector of the grid point
   *
   * @return level vector
   */
  const std::vector<level_t> getLevelVector() {
    if (level_vector.empty()) {
      for (size_t d = 0; d < point.getDimension(); d++) {
        level_vector.push_back(point.getLevel(d));
      }
    }

    return level_vector;
  }


  /***
   * Gets the grid point
   *
   * @return point
   */
  GridPoint& getPoint() {
    return point;
  }


  /***
   * Gets the sequential number of the grid point
   * @return sequential number
   */
  size_t getSeq() const {
    return seq;
  }

 private:
  GridPoint point;
  size_t seq;
  std::vector<level_t> level_vector;
};


/**
 * Abstract refinement class for sparse grids
 */
class AbstractRefinement {
 public:
  /**
  * Type of the identifier of the refinement atom (e.g. a grid point or a subspace)
  */
  typedef AbstractRefinement_refinement_key refinement_key_type;
  // typedef std::pair<size_t, size_t> refinement_key_type;


  /**
   * Type of functor value assigned to each refinement atom
   */
  typedef double refinement_value_type;  // refinement functor value


  /**
  * Pair for the refinement key and value
  */
  typedef std::pair<std::shared_ptr<refinement_key_type>, refinement_value_type>
  refinement_pair_type;

  typedef typename std::forward_list<AbstractRefinement::refinement_pair_type>
  refinement_list_type;


  /**
   * Comparison of the refinement_pair_type. This way the priority queue
   * has the elements with the smallest refinement_value_type on top
   */
  static bool compare_pairs(const refinement_pair_type& lhs,
                            const refinement_pair_type& rhs)  {
    return (lhs.second > rhs.second);
  }


  /**
  * Container for the collection of the refinement atoms and the corresponding
  * value
  */
  typedef std::vector<refinement_pair_type> refinement_container_type;


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
  virtual void free_refine(GridStorage& storage,
                           RefinementFunctor& functor,
                           std::vector<size_t>* addedPoints = nullptr) = 0;

  /**
   * Computes and returns the number of grid points, which can be refined.
   * This is the number of grid points that have at least one child missing.
   *
   * @param storage hashmap that stores the grid points
   * @return The number of grid points that can be refined
   */
  virtual size_t getNumberOfRefinablePoints(GridStorage& storage) = 0;


  /**
   * Refine one grid point along a single direction
   * @param storage hashmap that stores the grid points
   * @param point point to refine
   * @param d direction
   */
  virtual void refineGridpoint1D(GridStorage& storage, GridPoint& point,
                                 size_t d) = 0;


  /**
  * Refine one grid point along a single direction
  * @param storage hashmap that stores the grid points
  * @param seq sequential number of the grid point
  * @param d direction
  */
  virtual void refineGridpoint1D(GridStorage& storage, size_t seq, size_t d);


  /**
   * Check if the grid point is refinable
   * @param storage hashmap that stores the grid points
   * @param point grid point
   */
  bool isRefinable(GridStorage& storage, GridPoint& point);

  /**
   * Destructor
   */
  virtual ~AbstractRefinement() {
  }

  /**
   * Returns the index of the first occurrence of minimal element in array.
   * Used to find which entry is to be replaced next searching the maximum ones.
   *
   * @param array array with values
   * @param length length of array
   *
   * @return index of the first occurrence of minimal element in array
   */
  virtual size_t getIndexOfMin(RefinementFunctor::value_type* array,
                               size_t length);

 protected:
  /**
   * This method refines a grid point by generating the children in every dimension
   * of the grid and all their missing ancestors by calling create_gridpoint().
   *
   * @param storage hashmap that stores the gridpoints
   * @param refine_index The index in the hashmap of the point that should be refined
   */
  virtual void refineGridpoint(GridStorage& storage, size_t refine_index) = 0;

  /**
   * This method creates a new point on the grid. It checks if some parents or
   * children are needed in other dimensions.
   *
   * @param storage hashmap that stores the gridpoints
   * @param point The point that should be inserted
   */
  virtual void createGridpoint(GridStorage& storage, GridPoint& point) = 0;


  /**
   * Subroutine for grid point creation.
   *
   * @param storage hashmap that stores the gridpoints
   * @param point The point that should be inserted
   */
  virtual void createGridpointSubroutine(GridStorage& storage,
                                         GridPoint& point) {
    // For efficiency this function is defined the header file, this way it
    // be easily inlined by compiler.
    if (!storage.isContaining(point)) {
      // save old leaf value
      bool saveLeaf = point.isLeaf();
      point.setLeaf(false);
      createGridpoint(storage, point);
      // restore leaf value
      point.setLeaf(saveLeaf);
    } else {
      // set stored index to false
      storage.getPoint((storage.find(&point))->second).setLeaf(false);
    }
  }


  /**
   * Creates children grid points along single direction
   *
   * @param point The point that should be refined
   * @param d direction
   * @param storage hashmap that stores the gridpoints
   * @param source_index index value in the dimension d
   * @param source_level level value in the dimension d
   */
  virtual void createGridpoint1D(
      GridPoint& point,
    size_t d, GridStorage& storage,
    index_t& source_index, level_t& source_level);


  /**
   * Identifies the sparse grid refinement atoms (points or subspaces) with
   * the largest indicator values.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a RefinementFunctor specifying the refinement criteria
   * @param collection container with grid element identifiers (e.g. sequence number, grid point)
   *  and corresponding refinement values (usually empty)
   */
  virtual void collectRefinablePoints(
    GridStorage& storage,
    RefinementFunctor& functor,
    refinement_container_type& collection) = 0;


  /**
   * Refines the collection of points.
   *
   * @param storage hashmap that stores the grid points
   * @param functor a RefinementFunctor specifying the refinement criteria
   * @param collection container with grid element identifiers (e.g. sequence number, grid point)
   *  and corresponding refinement values
   */
  virtual void refineGridpointsCollection(
    GridStorage& storage,
    RefinementFunctor& functor,
    refinement_container_type& collection) = 0;


  /***
   * Gets a list of the elements with corresponding refinement indicators.
   *
   * The list usually contains only one element for refinement in all
   * direction and refinement indicators for individual direction for 1d
   * refinement
   *
   * @param storage hashmap that stores the grid points
   * @param iter grid_map iterator to the current grid point
   * @param functor refinement functor
   * @return
   */
  virtual refinement_list_type getIndicator(
    GridStorage& storage,
    const GridStorage::grid_map_iterator& iter,
    const RefinementFunctor& functor) const = 0;

  friend class
  // need to be a friend since it delegates the calls to
  // protected class methods
    RefinementDecorator;
};


}  // namespace base
}  // namespace sgpp

#endif /* ABSTRACTREFINEMENT_HPP */
