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
  std::vector<ScalarEvaluator> evaluators;
  std::vector<FloatArrayVector> basisCoefficients;
  std::vector<double> xValues;
  bool doesNeedParameter;

  void computeBasisCoefficients() {
    basisCoefficients = std::vector<FloatArrayVector>(xValues.size(), FloatArrayVector::zero());

    for (size_t i = 0; i < evaluators.size(); ++i) {
      auto coeff = evaluators[i].getBasisCoefficients();
      for (size_t j = 0; j < coeff.size(); ++j) {
        basisCoefficients[j].at(i) = coeff[j];
      }
    }
  }

 public:
  explicit ArrayEvaluator(bool doesNeedParameter,
                          ScalarEvaluator evaluatorPrototype = ScalarEvaluator())
      : evaluators(1, evaluatorPrototype),
        basisCoefficients(),
        xValues(),
        doesNeedParameter(doesNeedParameter) {}

  ArrayEvaluator(ArrayEvaluator<ScalarEvaluator> const &other)
      : evaluators(other.evaluators),
        basisCoefficients(other.basisCoefficients),
        xValues(other.xValues),
        doesNeedParameter(other.doesNeedParameter) {}

  ~ArrayEvaluator() {}

  virtual std::vector<FloatArrayVector> getBasisCoefficients() { return basisCoefficients; }

  virtual void setGridPoints(std::vector<double> const &xValues) {
    this->xValues = xValues;
    for (auto &eval : evaluators) {
      eval.setGridPoints(xValues);
    }

    computeBasisCoefficients();
  }

  virtual std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> cloneLinear() {
    return std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>(
        new ArrayEvaluator<ScalarEvaluator>(*this));
  }

  virtual bool needsOrderedPoints() { return evaluators[0].needsOrderedPoints(); }

  virtual bool needsParameter() { return doesNeedParameter; }

  virtual void setParameter(FloatArrayVector const &param) {
    if (evaluators.size() > param.size()) {
      evaluators.resize(param.size());
    } else {
      if (evaluators.size() == 0) {
        evaluators.emplace_back();
        evaluators.back().setGridPoints(xValues);
      }
      while (evaluators.size() < param.size()) {
        evaluators.push_back(evaluators.back());
      }
    }
    for (size_t i = 0; i < param.size(); ++i) {
      evaluators[i].setParameter(param[i]);
    }

    computeBasisCoefficients();
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_ARRAYEVALUATOR_HPP_ */
