// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDGENERATOR_HPP
#define GRIDGENERATOR_HPP

#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <set>
#include <unordered_set>
#include <vector>

namespace sgpp {
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
   * Creates a regular sparse grid for a certain level @f$ n @f$, i.e., @f$ V_n^{(1)} =
   *\bigoplus_{|\vec{l}|_1 \leq n+d-1} W_{\vec{l}}@f$.
   *
   * @param level Grid level
   */
  virtual void regular(size_t level) = 0;

  /**
   * Creates a regular sparse grid for a certain level @f$ n @f$, i.e., @f$ V_n^{(1)} =
   *\bigoplus_{|\vec{l}|_1 \leq n+d-1} W_{\vec{l}}@f$.
   * If the used grid doesn't support the parameter t, t = 0 is used instead.
   *
   * @param level Grid level
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces
   */
  virtual void regular(size_t level, double T) {
    throw sgpp::base::not_implemented_exception("Parameter T not implemented for this grid type!");
  }

  /**
   * Creates a regular sparse grid for a certain level @f$ n @f$, i.e., @f$ V_n^{(1)} =
   *\bigoplus_{|\vec{l}|_1 \leq n+d-1} W_{\vec{l}}@f$.
   * If the used grid doesn't support the parameter t, t = 0 is used instead.
   * This grid generator allows the creation of regular grids that only contain
   * some interaction terms.
   * @param level Grid level
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces
   * @param terms determines the included interaction terms.
   */
  virtual void regularInter(size_t level, const std::set<std::set<size_t>>& terms, double T) {
    throw sgpp::base::not_implemented_exception(
        "Interaction-Term aware sparse grids not implemented for this grid type!");
  }

  /**
   * Creates a sparse grid with fully connected cliques
   *
   * @param level Grid level
   * @param clique_size clique size
   */
  virtual void cliques(size_t level, size_t clique_size) = 0;

  /**
   * Creates a sparse grid with fully connected cliques
   *
   * @param level Grid level
   * @param clique_size clique size
   * @param T modifier for subgrid selection, T = 0 implies standard sparse grid.
   *        For further information see Griebel and Knapek's paper
   *        optimized tensor-product approximation spaces
   */
  virtual void cliques(size_t level, size_t clique_size, double T) {
    throw sgpp::base::not_implemented_exception("Parameter T not implemented for this grid type!");
  }

  /**
   * Creates a full grid for a certain level @f$ n @f$, i.e., @f$ V_n = \bigoplus_{|\vec{l}|_\infty
   *\leq n} W_{\vec{l}}@f$.
   *
   * @param level Grid level
   */
  virtual void full(size_t level) = 0;

  /**
   * Creates an anisotropicFull Grid for certain level vector, for example [2,3,1] results
   * in a 3x7x1 grid
   */

  virtual void anisotropicFull(std::vector<size_t> dimlevels) {}

  /**
   * Creates a grid which doesn't contain the fullgrids with li<l_user, for any li level_t
   * */

  virtual void truncated(size_t level, size_t l_user) {}

  /**
   * Refines a grid according to the settings of the RefinementFunctor func.
   *
   * @param func reference to refinement functor
   * @param addedPoints pointer to vector to add newly created grid points to
   */
  virtual void refine(RefinementFunctor& func, std::vector<size_t>* addedPoints = nullptr) = 0;

  /**
   * Refines a grid according to the settings of the RefinementFunctor func.
   * Does not create any interactions, that are not in the list of allowed interactions.
   * @details Refines the grid, but only adds interactions that are contained in
   * the set interactions, i.e. only desired interactions.
   * Each desired interaction is encoded as a vector which contains all desired interactions.
   * For example, if we want to include grid points that model an
   * interaction between the first and the second predictor, we would
   * include the vector [1,2] in interactions.
   * @param func pointer to refinement functor
   * @param interactions allowed interactions
   */
  virtual void refineInter(RefinementFunctor& func,
                           const std::set<std::set<size_t>>& interactions) {
    throw sgpp::base::not_implemented_exception(
        "Interaction-Term refinement not implemented for this grid type!");
  }
  /**
   * Coarsens a  grid according to the settings of the CoarseningFunctor func.
   *
   * @param func pointer to coarsening functor
   * @param alpha Pointer to DataVector containing the grid's coefficients
   * @param removedSeq pointer to vector to append the seq numbers of coarsened grid points to.
   */
  virtual void coarsen(CoarseningFunctor& func, DataVector& alpha,
                       std::vector<size_t>* removedSeq) = 0;

  /**
   * Coarsens a  grid according to the settings of the CoarseningFunctor func.
   * Only numFirstOnly first grid points are checked for coarsening.
   *
   * @param func pointer to coarsening functor
   * @param alpha Pointer to DataVector containing the grid's coefficients
   * @param numFirstOnly max. number grid points to be coarsened
   * @param removedSeq pointer to vector to append the seq numbers of coarsened grid points to.
   */
  virtual void coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha, size_t numFirstOnl,
                                 std::vector<size_t>* removedSeqy) = 0;

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
  virtual void refineMaxLevel(RefinementFunctor& func, size_t maxLevel) = 0;

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
}  // namespace sgpp

#endif /* GRIDGENERATOR_HPP */
