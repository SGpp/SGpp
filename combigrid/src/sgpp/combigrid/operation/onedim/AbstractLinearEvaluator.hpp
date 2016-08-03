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
 * Linear means that the evaluator is linear in the function it interpolates, not necessarily (in
 * the interpolation case) in the evaluation point x.
 */
template <typename V>
class AbstractLinearEvaluator : public AbstractEvaluator<V> {
 public:
  virtual ~AbstractLinearEvaluator() {}

  virtual std::vector<V> getBasisCoefficients() = 0;

  virtual void setGridPoints(std::vector<double> const &xValues) = 0;
  virtual std::shared_ptr<AbstractLinearEvaluator<V>> cloneLinear() = 0;
  virtual std::shared_ptr<AbstractEvaluator<V>> clone() { return cloneLinear(); }
  virtual bool needsOrderedPoints() = 0;
  virtual bool needsParameter() = 0;
  virtual void setParameter(V const &param) = 0;

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
