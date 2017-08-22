// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class takes a 1D linear evaluator operating on scalars (which is given through the template
 * type ScalarEvaluator) and uses multiple instances of it for doing multi-evaluation.
 * It is optimized to perform operations common to all of the evaluators only once.
 */
template <typename ScalarEvaluator>
class ArrayEvaluator : public AbstractLinearEvaluator<FloatArrayVector> {
  ScalarEvaluator evaluator;
  FloatArrayVector params;
  std::vector<FloatArrayVector> basisCoefficients;
  bool coefficientsComputed = false;
  std::vector<double> xValues;
  bool doesNeedParameter;

  void computeBasisCoefficients() {
    basisCoefficients = std::vector<FloatArrayVector>(xValues.size(), FloatArrayVector::zero());

    if (params.size() == 0) {
      auto coeff = evaluator.getBasisCoefficients();
      for (size_t j = 0; j < coeff.size(); ++j) {
        basisCoefficients[j].at(0) = coeff[j];
      }
    } else {
      for (size_t i = 0; i < params.size(); ++i) {
        evaluator.setParameter(params[i]);
        auto coeff = evaluator.getBasisCoefficients();
        for (size_t j = 0; j < coeff.size(); ++j) {
          basisCoefficients[j].at(i) = coeff[j];
        }
      }
    }

    coefficientsComputed = true;
  }

 public:
  explicit ArrayEvaluator(bool doesNeedParameter,
                          ScalarEvaluator evaluatorPrototype = ScalarEvaluator())
      : evaluator(evaluatorPrototype),
        params(),
        basisCoefficients(),
        xValues(),
        doesNeedParameter(doesNeedParameter) {}

  explicit ArrayEvaluator(bool doesNeedParameter, ScalarEvaluator *evaluatorPrototype)
      : evaluator(*evaluatorPrototype),
        params(),
        basisCoefficients(),
        xValues(),
        doesNeedParameter(doesNeedParameter) {}

  ArrayEvaluator(ArrayEvaluator<ScalarEvaluator> const &other)
      : evaluator(other.evaluator),
        params(other.params),
        basisCoefficients(other.basisCoefficients),
        coefficientsComputed(other.coefficientsComputed),
        xValues(other.xValues),
        doesNeedParameter(other.doesNeedParameter) {}

  ~ArrayEvaluator() {}

  virtual std::vector<FloatArrayVector> getBasisCoefficients() {
    if (!coefficientsComputed) {
      computeBasisCoefficients();
    }
    return basisCoefficients;
  }

  virtual void setGridPoints(std::vector<double> const &xValues) {
    this->xValues = xValues;
    evaluator.setGridPoints(xValues);

    coefficientsComputed = false;
  }

  virtual std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> cloneLinear() {
    return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>(
        new ArrayEvaluator<ScalarEvaluator>(*this));
  }

  virtual bool needsOrderedPoints() { return evaluator.needsOrderedPoints(); }

  virtual bool needsParameter() { return doesNeedParameter; }

  virtual void setParameter(FloatArrayVector const &param) {
    this->params = param;
    coefficientsComputed = false;
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_ */
