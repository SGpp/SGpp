// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTLINEAREVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTLINEAREVALUATOR_HPP_

#include <sgpp/combigrid/operation/onedim/AbstractEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * See also AbstractEvaluator.
 * Linear means that the evaluator is linear in the function it interpolates, not necessarily (in
 * the interpolation case) in the evaluation point x.
 * Especially, this means, that linear coefficients for the function values at each grid points can
 * be calculated, which are used by FullGridTensorEvaluator.
 */
template <typename V>
class AbstractLinearEvaluator : public AbstractEvaluator<V> {
 public:
  virtual ~AbstractLinearEvaluator() {}

  /**
   * This method is specific to AbstractLinearEvaluator and has to be implemented by all subclasses.
   * It should return the coefficients used for the linear combination of function values.
   * Example: If you want to do a linear interpolation and have the grid points 0.0 and 0.5 and an
   * evaluation point (the parameter) of 0.2, you would return the weights (0.6, 0.4), because the
   * interpolated value is 0.6 * f(0) + 0.4 * f(0.5).
   *
   * In the case of multi-evaluation, this returns several coefficients for each grid point.
   * Continuing the example above, if you would have evaluation points 0.2 and 0.4, the basis
   * coefficients returned should be ((0.6, 0.2), (0.4, 0.8)).
   */
  virtual std::vector<V> getBasisCoefficients() = 0;

  virtual void setGridPoints(std::vector<double> const &xValues) = 0;

  /**
   * Clones this object and returns it as a shared pointer to AbstractLinearEvaluator<V>.
   */
  virtual std::shared_ptr<AbstractLinearEvaluator<V>> cloneLinear() = 0;
  virtual std::shared_ptr<AbstractEvaluator<V>> clone() { return cloneLinear(); }
  virtual bool needsOrderedPoints() = 0;
  virtual bool needsParameter() = 0;
  virtual void setParameter(V const &param) = 0;

  /**
   * AbstractLinearEvaluator provides a standard implementation of this method based on
   * getBasisCoefficients().
   */
  virtual V eval(std::vector<double> const &functionValues) {
    auto basisCoefficients = getBasisCoefficients();

    V sum = V::zero();
    for (size_t i = 0; i < functionValues.size(); ++i) {
      V v = basisCoefficients[i];
      v.scalarMult(functionValues[i]);
      sum.add(v);
    }

    return sum;
  }

  /**
   * AbstractLinearEvaluator provides a standard implementation of this method based on
   * getBasisCoefficients().
   */
  virtual V eval(std::vector<V> const &functionValues) {
    auto basisCoefficients = getBasisCoefficients();

    V sum = V::zero();

    for (size_t i = 0; i < functionValues.size(); ++i) {
      V v = basisCoefficients[i];
      v.componentwiseMult(functionValues[i]);
      sum.add(v);
    }

    return sum;
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTLINEAREVALUATOR_HPP_ */
