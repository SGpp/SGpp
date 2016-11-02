// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_

#include <sgpp/globaldef.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Abstract base class for one-dimensional numerical evaluation methods. Currently, only its derived
 * class AbstractLinearEvaluator is used.
 *
 * The template parameter is used for switching between single-evaluation and multi-evaluation, see
 * FloatArrayVector.
 */
template <typename V>
class AbstractEvaluator {
 public:
  virtual ~AbstractEvaluator() {}

  /**
   * Sets the grid points for evaluation. The evaluation method should allow for arbitrary grid
   * points to achieve maximal modularity.
   */
  virtual void setGridPoints(std::vector<double> const &xValues) = 0;

  /**
   * Clones the AbstractEvaluator object.
   */
  virtual std::shared_ptr<AbstractEvaluator<V>> clone() = 0;

  /**
   * @return true iff the grid points have to be provided in ascending order (for example for linear
   * interpolation).
   */
  virtual bool needsOrderedPoints() = 0;

  /**
   * @return true iff this evaluation method needs a parameter (e.g. true for interpolation, false
   * for quadrature).
   */
  virtual bool needsParameter() = 0;

  /**
   * Via this method, the parameter can be set.
   */
  virtual void setParameter(V const &param) = 0;

  /**
   * Evaluates the numerical method on the already set grid points.
   * @param functionValues function values at the grid points, the order must match the one of the
   * grid points.
   */
  virtual V eval(std::vector<double> const &functionValues) = 0;

  /**
   * (Currently not used). Via the template type V, multiple arrays of function values for multiple
   * functions might be specified and a result for each of these functions should be returned.
   * This could be used in an alternative implementation to FullGridTensorEvaluator that would also
   * be able to handle non-linear one-dimensional methods.
   */
  virtual V eval(std::vector<V> const &functionValues) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTEVALUATOR_HPP_ */
