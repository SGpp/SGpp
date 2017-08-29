// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTPRECALCULATEDWEIGHTSEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTPRECALCULATEDWEIGHTSEVALUATOR_HPP_

#include <sgpp/combigrid/operation/onedim/AbstractEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * See also AbstractEvaluator.
 * Precalculated coefficients  mean that in contrary to AbstractLinearEvaluator not simply the
 * function values are weighted by weights depending on the evaluation point but that more advanced
 * basis functions (for example B-Splines) are used for which the interpolation weights need to
 * be precalculated via the method ToDo(rehmemk) precalculateBsplineweights
 */
template <typename V>
class AbstractPrecalculatedWeightsEvaluator : public AbstractEvaluator<V> {
 public:
  virtual ~AbstractPrecalculatedWeightsEvaluator() {}

  /**
   * This method is specific to AbstractPrecalculatedWeightsEvaluator and has to be implemented
   * by all subclasses.
   * It should precalculate the weights used for the linear combination of basis function values.
   *
   * In the case of multi-evaluation, this returns several coefficients for each grid point.
   */
  virtual std::vector<V> getBasisWeights() = 0;

  //??? Do we need this for B-Splines. Should be covered by the prepare method
  virtual void setGridPoints(std::vector<double> const &xValues) = 0;

  /**
   * Clones this object and returns it as a shared pointer to
   * AbstractPrecalculatedCoefficientsEvaluator<V>.
   */
  virtual std::shared_ptr<AbstractPrecalculatedWeightsEvaluator<V>> clonePrecalculatedWeights() = 0;
  virtual std::shared_ptr<AbstractEvaluator<V>> clone() { return clonePrecalculatedWeights(); }

  // ??? Does B-Spline evaluation need ordered points?
  virtual bool needsOrderedPoints() = 0;
  virtual bool needsParameter() = 0;
  virtual void setParameter(V const &param) = 0;

  /**
   * AbstractPrecalculatedWeightsEvaluator provides a standard implementation of this method
   * based on getBasisWeights().
   */

  // ToDo(rehmemk) implement B-Splines and method to retrieve BasisfunctionValues
  virtual V eval(std::vector<double> const &BasisfunctionValues) {
    auto basisWeights = getBasisWeights();

    V sum = V::zero();
    for (size_t i = 0; i < BasisfunctionValues.size(); ++i) {
      V v = basisWeights[i];
      v.scalarMult(BasisfunctionValues[i]);
      sum.add(v);
    }

    return sum;
  }

  /**
   * AbstractLinearEvaluator provides a standard implementation of this method based on
   * getBasisCoefficients().
   */
  virtual V eval(std::vector<V> const &BasisfunctionValues) {
    auto basisWeights = getBasisWeights();

    V sum = V::zero();

    for (size_t i = 0; i < BasisfunctionValues.size(); ++i) {
      V v = basisWeights[i];
      v.componentwiseMult(BasisfunctionValues[i]);
      sum.add(v);
    }

    return sum;
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ABSTRACTPRECALCULATEDWEIGHTSEVALUATOR_HPP_ \
          */
