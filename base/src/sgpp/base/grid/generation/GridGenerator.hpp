// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDGENERATOR_HPP
#define GRIDGENERATOR_HPP

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Abstract class that defines the interfaces for the different grid's GridGenerators
 */
class GridGenerator {
 public:
  /**
   * Constructor
   */
  GridGenerator() {}

  /**
   * Destructor
   */
  virtual ~GridGenerator() {}

  /**
   * Creates a regular sparse grid for a certain level @f$ n @f$, i.e., @f$ V_n^{(1)} = \bigoplus_{|\vec{l}|_1 \leq n+d-1} W_{\vec{l}}@f$.
   *
   * @param level Grid level
   */
  virtual void regular(size_t level) = 0;


  /**
   * Creates a sparse grid with fully connected cliques
   *
   * @param level Grid level
   * @param clique_size clique size
   */
  virtual void cliques(size_t level, size_t clique_size) = 0;

  /**
   * Creates a full grid for a certain level @f$ n @f$, i.e., @f$ V_n = \bigoplus_{|\vec{l}|_\infty \leq n} W_{\vec{l}}@f$.
   *
   * @param level Grid level
   */
  virtual void full(size_t level) = 0;

  /**
   * Creates a grid which doesn't contain the fullgrids with li<l_user, for any li level_t
   * */
  virtual void truncated(size_t level, size_t l_user) {}

  /**
   * Refines a grid according to the settings of the RefinementFunctor func.
   *
   * @param func pointer to refinement functor
   */
  virtual void refine(RefinementFunctor* func) = 0;

  /**
   * Coarsens a  grid according to the settings of the CoarseningFunctor func.
   *
   * @param func pointer to coarsening functor
   * @param alpha Pointer to DataVector containing the grid's coefficients
   */
  virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) = 0;

  /**
   * Coarsens a  grid according to the settings of the CoarseningFunctor func.
   * Only numFirstOnly first grid points are checked for coarsening.
   *
   * @param func pointer to coarsening functor
   * @param alpha Pointer to DataVector containing the grid's coefficients
   * @param numFirstOnly max. number grid points to be coarsened
   */
  virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha,
                                 size_t numFirstOnly) = 0;

  /**
   * Returns the number of points on the grid that can be refined in the next iteration
   *
   * @return the number of points on the grid that can be refined
   */
  virtual size_t getNumberOfRefinablePoints() = 0;

  /**
   * Returns the number of points on the grid that can be removed in the next iteration
   *
   * @return the number of points on the grid that can be removed
   */
  virtual size_t getNumberOfRemovablePoints() = 0;

  /**
   * Refines a grid according to the settings of the RefinementFunctor func.
   * additionally a maximum level for refinement is taken into account
   *
   * @param func pointer to refinement functor
   * @param maxLevel no points on higher levels than maxLevel will be created
   */
  virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) = 0;

  /**
   * Returns the number of points on the grid that can be refined in the next iteration
   * additionally a maximum level for refinement is taken into account
   *
   * @param maxLevel no points on higher levels than maxLevel will be created
   *
   * @return the number of points on the grid that can be refined
   */
  virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) = 0;
};

}  // namespace base
}  // namespace SGPP

#endif /* GRIDGENERATOR_HPP */
