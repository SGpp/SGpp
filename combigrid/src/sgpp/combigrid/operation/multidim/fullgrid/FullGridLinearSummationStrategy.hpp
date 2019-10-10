// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <iostream>
#include <vector>

namespace sgpp {
namespace combigrid {

template <typename V>
class FullGridLinearSummationStrategy : public AbstractFullGridSummationStrategy<V> {
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
  FullGridLinearSummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, evaluatorPrototypes, pointHierarchies) {}

  ~FullGridLinearSummationStrategy() override {}

  /**
   * Evaluates the function given through the storage for a certain level-multi-index (see class
   * description).
   * Summation of the form \f$\sum_i \alpha_i basis_i(param) \f$
   * This is used for interpolation and quadratures
   */
  V eval(MultiIndex const &level) override {
    CGLOG("FullGridTensorEvaluator::eval(): start");
    size_t numDimensions = this->evaluators.size();
    size_t lastDim = numDimensions - 1;
    MultiIndex multiBounds(numDimensions);
    std::vector<bool> orderingConfiguration(numDimensions);

    size_t paramIndex = 0;

    // the basis coefficients for this level are stored
    // the bounds for traversal are initialized given the number of points in each direction
    // if not already stored, the evaluators for the given level are cloned and their parameter, if
    // needed, is set
    // the ordering configuration is created - it stores a boolean for each dimension expressing
    // whether the corresponding evaluator needs sorted (ascending) points

    // init evaluators and basis values, init multiBounds and orderingConfiguration
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
    // init partial products
    this->partialProducts[0] = V::one();
    for (size_t d = 1; d < numDimensions; ++d) {
      V value = this->partialProducts[d - 1];
      value.componentwiseMult(this->basisValues[d - 1][0]);
      this->partialProducts[d] = value;
    }

    CGLOG("FullGridTensorEvaluator::eval(): create storage iterator");

    // start iteration
    MultiIndexIterator it(multiBounds);

    auto funcIter = this->storage->getGuidedIterator(level, it, orderingConfiguration);
    V sum = V::zero();
    CGLOG("FullGridTensorEvaluator::eval(): start loop");

    if (!funcIter->isValid()) {  // should not happen
      return sum;
    }

    while (true) {
      CGLOG("FullGridTensorEvaluator::eval(): in loop");
      // get function value and partial product and multiply them together with the last basis
      // coefficient, then add the resulting value to the total sum
      double value = funcIter->value();
      //      std::cout << std::defaultfloat << value << " ";
      V vec = this->partialProducts[lastDim];
      vec.componentwiseMult(this->basisValues[lastDim][it.indexAt(lastDim)]);
      vec.scalarMult(value);
      sum.add(vec);

      // increment iterator
      int h = funcIter->moveToNext();

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
    //    std::cout << "\n";
    return sum;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
