// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

enum class FullGridSummationStrategyType {
  LINEAR,
  QUADRATIC,
  VARIANCE,
  TENSORVARIANCE,
  FULLSUBSPACEDPCE,
  ONEDSUBSPACEPCE
};

template <typename V>
class AbstractFullGridSummationStrategy {
 protected:
  std::shared_ptr<AbstractCombigridStorage> storage;
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;

  // partialProducts[i] stores the product of the first i basis values (corresponding to the current
  // multi-index) , i. e. partialProducts[0] = 1
  // partialProducts has Size numDimensions, since the product partialProducts[numDimensions] is
  // only used once and does not have to be stored
  std::vector<V> partialProducts;

  /**
   * For each dimension, this contains a vector of weights which are used as coefficients for
   * linearly combining the function values at different grid points.
   */
  std::vector<std::vector<V>> basisValues;

  // one per dimension and level
  std::vector<std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>>> evaluators;
  // parameters (empty when doing quadrature)
  std::vector<V> parameters;

 public:
  /**
   * Constructor.
   *
   * @param storage Storage that stores and provides the function values for each grid point.
   * @param evaluatorPrototypes prototype objects for the evaluators that are cloned to get an
   * evaluator for each dimension and each level.
   * @param pointHierarchies PointHierarchy objects for each dimension providing the points for each
   * level and information about their ordering.
   */
  AbstractFullGridSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : storage(storage),
        evaluatorPrototypes(evaluatorPrototypes),
        pointHierarchies(pointHierarchies),
        partialProducts(evaluatorPrototypes.size()),
        basisValues(evaluatorPrototypes.size()),
        evaluators(evaluatorPrototypes.size()),
        parameters(evaluatorPrototypes.size()) {
    // TODO(holzmudd): check for dimension equality
  }

  virtual ~AbstractFullGridSummationStrategy() {}

  virtual V eval(MultiIndex const &level) = 0;

  /**
   * Sets the parameters for the evaluators. Each dimension in which the evaluator does not need a
   * parameter is skipped.
   * So if only the evaluators at dimensions 1 and 3 need a parameter, params.size() should be 2 (or
   * at least 2)
   */
  void setParameters(std::vector<V> const &params) {
    size_t numDimensions = this->evaluatorPrototypes.size();

    parameters = params;

    size_t paramIndex = 0;

    // we can't just set the parameters to the prototypes because the prototypes might be identical
    // (the pointer to one prototype might be duplicated)
    for (size_t d = 0; d < numDimensions; ++d) {
      auto &prototype = this->evaluatorPrototypes[d];

      if (prototype->needsParameter()) {
        if (paramIndex >= params.size()) {
          throw std::runtime_error(
              "AbstractFullGridSummationstrategy::setParameters(): parameter dimensionality is too "
              "low.");
        }
        // prototype->setParameter(params[paramIndex]); <- this is useless, see above
        for (auto &eval : evaluators[d]) {
          eval->setParameter(params[paramIndex]);
        }

        ++paramIndex;
      }
    }

    // can occur for quadrature, since a "default parameter"
    // may be passed
    /*if (paramIndex < params.size()) {
      throw std::runtime_error(
          "AbstractFullGridLinearEvaluator::setParameters(): parameter dimensionality is too "
          "high.");
    } */
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
