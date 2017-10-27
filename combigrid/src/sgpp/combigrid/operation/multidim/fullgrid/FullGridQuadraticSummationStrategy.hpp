// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <iostream>
#include <vector>
#include "AbstractFullGridEvaluationStrategy.hpp"

namespace sgpp {
namespace combigrid {

template <typename V>
class FullGridQuadraticSummationStrategy : public AbstractFullGridSummationStrategy<V> {
 protected:
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
  FullGridQuadraticSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, evaluatorPrototypes, pointHierarchies) {}

  FullGridQuadraticSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      GridFunction gridFunction)
      : AbstractFullGridSummationStrategy<V>(storage, evaluatorPrototypes, pointHierarchies,
                                             gridFunction) {}

  ~FullGridQuadraticSummationStrategy() {}

  /**
   * Evaluates the function given through the storage for a certain level-multi-index (see class
   * description).
   * Summation of the form \sum_i \alpha_i \sum_j \alpha_j basis_ij(x)
   * This is used for variance calculations where basis_ij = \int basis_i(x) * basis_j(x) dx
   */
  V eval(MultiIndex const &level) {
    if (this->strategy == Strategy::grid_based) {
      if (!this->precomputedLevels->containsIndex(level)) {
        this->addResults(level, this->gridFunction(this->getTensorGrid2(level)));
        this->precomputedLevels->set(level, 1);
      }
    }
    CGLOG("FullGridTensorEvaluator::eval(): start");
    size_t numDimensions = this->evaluators.size();
    size_t lastDim = numDimensions - 1;
    MultiIndex multiBounds(numDimensions);
    std::vector<bool> orderingConfiguration(numDimensions);

    size_t paramIndex = 0;

    // *the basis coefficients for this level are stored
    // *the bounds for traversal are initialized given the number of points in each direction
    // *if not already stored, the evaluators for the given level are cloned and their parameter, if
    //  needed, is set
    // *the ordering configuration is created - it stores a boolean for each dimension expressing
    //  whether the corresponding evaluator needs sorted (ascending) points

    // initialize evaluators and basis values, initialize multiBounds and orderingConfiguration
    for (size_t d = 0; d < numDimensions; ++d) {
      size_t currentLevel = level[d];
      auto &currentEvaluators = this->evaluators[d];

      bool needsParam = this->evaluatorPrototypes[d]->needsParameter();
      bool needsOrdered = this->evaluatorPrototypes[d]->needsOrderedPoints();

      for (size_t l = currentEvaluators.size(); l <= currentLevel; ++l) {
        auto eval = this->evaluatorPrototypes[d]->cloneLinear();
        eval->setGridPoints(this->pointHierarchies[d]->getPoints(l, needsOrdered));
        eval->setLevel(l);
        if (needsParam) {
          eval->setParameter(this->parameters[paramIndex]);
        }
        currentEvaluators.push_back(eval);
      }
      this->basisValues[d] = currentEvaluators[currentLevel]->getBasisValues();
      multiBounds[d] = this->pointHierarchies[d]->getNumPoints(currentLevel);
      orderingConfiguration[d] = needsOrdered;

      if (needsParam) {
        ++paramIndex;
      }
    }

    // for efficient computation, the products over the first i evaluator coefficients are stored
    // for all i up to n-1.
    // This way, we only have to multiply them with the values for the changing indices at each
    // iteration step.

    // initialize partial products
    this->partialProducts[0] = V::one();
    for (size_t d = 1; d < numDimensions; ++d) {
      V value = this->partialProducts[d - 1];
      value.componentwiseMult(this->basisValues[d - 1][0]);
      this->partialProducts[d] = value;
    }

    CGLOG("FullGridTensorEvaluator::eval(): create storage iterator");

    // start iteration
    MultiIndexIterator it(multiBounds);

    auto funcIter_j = this->storage->getGuidedIterator(level, it, orderingConfiguration);
    V sum = V::zero();

    CGLOG("FullGridTensorEvaluator::eval(): start loop");

    if (!funcIter_j->isValid()) {  // should not happen
      return sum;
    }

    // ToDo (rehmemk) Finish this!s

    // get partial product and multiply them together with the last basis
    // coefficient, multiply them with the i-coefficient and add the resulting value to the inner
    // sum. Then get the j-coefficient, multiply it by the inner sum and add this to (the outer) sum
    while (true) {
      CGLOG("FullGridTensorEvaluator::eval(): in loop");

      double value_j = funcIter_j->value();
      V vec = this->partialProducts[lastDim];
      vec.componentwiseMult(this->basisValues[lastDim][it.indexAt(lastDim)]);
      vec.scalarMult(value_j);
      sum.add(vec);

      // increment iterator
      int h = funcIter_j->moveToNext();

      CGLOG("FullGridTensorEvaluator::eval(): moveToNext() == " << h);

      // check if not only the last dimension index changed
      if (h != 0) {
        if (h < 0) {
          break;  // all indices have been traversed, stop iteration and return sum
        } else {
          // more than the last index have changed, thus update partialProducts
          for (size_t d = lastDim - h; d < lastDim; ++d) {
            auto pp = this->partialProducts[d];  // TODO(holzmudd): could probably be optimized...
            pp.componentwiseMult(this->basisValues[d][it.indexAt(d)]);
            this->partialProducts[d + 1] = pp;
          }
        }
      }
    }

    it.reset();
    auto funcIter_i = this->storage->getGuidedIterator(level, it, orderingConfiguration);
    // Multiply sum[i] by alpha_i
    for (size_t i = 0; i < sum.size(); i++) {
      double value_i = funcIter_i->value();
      sum[i] *= value_i;
      funcIter_i->moveToNext();
    }
    return sum;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
