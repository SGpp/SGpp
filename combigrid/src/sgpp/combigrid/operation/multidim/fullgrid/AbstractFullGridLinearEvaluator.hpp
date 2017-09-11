// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <iostream>
#include <memory>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FunctionValuesCoefficientsStorage.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/SLECoefficientsStorage.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

template <typename V>
class AbstractFullGridLinearEvaluator : public AbstractFullGridEvaluator<V> {
 protected:
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
  MultiIndex multiBounds;

  // one per dimension
  std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes;
  // one per dimension and level
  std::vector<std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>>> evaluators;
  // parameters (empty when doing quadrature)
  std::vector<V> parameters;
  // recompute basis values flag
  bool updateBasisValues;

  // storage for coefficients of linear combination
  std::shared_ptr<AbstractBasisCoefficientsStorage> coefficientsStorage;

  // store ordering configuration of one-dimensional evaluators
  std::vector<bool> orderingConfiguration;

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
  AbstractFullGridLinearEvaluator(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridEvaluator<V>(storage, pointHierarchies),
        partialProducts(evaluatorPrototypes.size()),
        basisValues(evaluatorPrototypes.size()),
        multiBounds(evaluatorPrototypes.size()),
        evaluatorPrototypes(evaluatorPrototypes),
        evaluators(evaluatorPrototypes.size()),
        parameters(evaluatorPrototypes.size()),
        updateBasisValues(false),
        orderingConfiguration(evaluatorPrototypes.size()) {
    // TODO(holzmudd): check for dimension equality

    // init orderingConfiguration
    for (size_t d = 0; d < evaluators.size(); ++d) {
      orderingConfiguration[d] = evaluatorPrototypes[d]->needsOrderedPoints();
    }

    initializeCoefficientStorage();
  }

  virtual ~AbstractFullGridLinearEvaluator() {}

  void initializeCoefficientStorage() {
    size_t d = 0;
    while (d < evaluatorPrototypes.size()) {
      auto &currentEvaluator = evaluatorPrototypes[d];
      if (currentEvaluator->getBasisCoefficientComputationType() ==
          BasisCoefficientsComputationType::SLE) {
        coefficientsStorage = std::make_shared<SLECoefficientsStorage>();
      }
      d++;
    }

    coefficientsStorage = std::make_shared<FunctionValuesCoefficientsStorage>();
  }

  /**
   * Evaluates the function given through the storage for a certain level-multi-index (see class
   * description).
   */
  virtual V eval(MultiIndex const &level) {
    CGLOG("FullGridTensorEvaluator::eval(): start");
    size_t numDimensions = evaluators.size();
    size_t lastDim = numDimensions - 1;

    // evaluate the basis functions of the subspace specified by the level variable
    // at the current evaluation points
    computeBasisValues(level);

    // for efficient computation, the products over the first i evaluator coefficients are stored
    // for all i up to n-1.
    // This way, we only have to multiply them with the values for the changing indices at each
    // iteration step.
    // init partial products
    partialProducts[0] = V::one();
    for (size_t d = 1; d < numDimensions; ++d) {
      V value = partialProducts[d - 1];
      value.componentwiseMult(basisValues[d - 1][0]);
      partialProducts[d] = value;
    }

    CGLOG("FullGridTensorEvaluator::eval(): create storage iterator");

    // prepare coefficients
    std::shared_ptr<std::vector<double>> coefficients = coefficientsStorage->getCoefficients(
        level, this->storage, multiBounds, orderingConfiguration);

    // start iteration
    MultiIndexIterator it(multiBounds);
    auto funcIter = this->storage->getGuidedIterator(level, it, orderingConfiguration);
    V sum = V::zero();
    if (!funcIter->isValid()) {  // should not happen
      return sum;
    }

    CGLOG("FullGridTensorEvaluator::eval(): start loop");
    size_t i = 0;
    while (true) {
      CGLOG("FullGridTensorEvaluator::eval(): in loop");
      // get function value and partial product and multiply them together with the last basis
      // coefficient, then add the resulting value to the total sum
      double coefficient = (*coefficients)[i++];
      V vec = partialProducts[lastDim];
      vec.componentwiseMult(basisValues[lastDim][it.indexAt(lastDim)]);
      vec.scalarMult(coefficient);
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
            auto pp = partialProducts[d];  // TODO(holzmudd): could probably be optimized...
            pp.componentwiseMult(basisValues[d][it.indexAt(d)]);
            partialProducts[d + 1] = pp;
          }
        }
      }
    }

    return sum;
  }

  /**
   * Sets the parameters for the evaluators. Each dimension in which the evaluator does not need a
   * parameter is skipped.
   * So if only the evaluators at dimensions 1 and 3 need a parameter, params.size() should be 2 (or
   * at least 2)
   */
  void setParameters(std::vector<V> const &params) {
    size_t numDimensions = evaluatorPrototypes.size();

    parameters = params;
    updateBasisValues = true;

    size_t paramIndex = 0;

    // we can't just set the parameters to the prototypes because the prototypes might be identical
    // (the pointer to one prototype might be duplicated)
    for (size_t d = 0; d < numDimensions; ++d) {
      auto &prototype = evaluatorPrototypes[d];

      if (prototype->needsParameter()) {
        if (paramIndex >= params.size()) {
          throw std::runtime_error(
              "AbstractFullGridLinearEvaluator::setParameters(): parameter dimensionality is too "
              "low ");
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

  void computeBasisValues(MultiIndex const &level) {
    size_t numDimensions = evaluators.size();

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
      auto &currentEvaluators = evaluators[d];

      bool needsParam = evaluatorPrototypes[d]->needsParameter();
      for (size_t l = currentEvaluators.size(); l <= currentLevel; ++l) {
        auto eval = evaluatorPrototypes[d]->cloneLinear();
        eval->setLevel(l);

        // ToDo (rehmemk) Add evaluated B-Spline values to GridPoints in method setGridPoints in
        // LinearInterpolationEvaluator.
        BasisCoefficientsComputationType new_basistype = BasisCoefficientsComputationType::SLE;
        eval->setBasisCoefficientComputationType(new_basistype);
        eval->setGridPoints(this->pointHierarchies[d]->getPoints(l, orderingConfiguration[d]));

        // ToDo (rehmemk) Add SLE solving to setFunctionValuesAtGridPoints in
        // LinearInterpolationEvaluator.cpp and call it here so that as 'FunctionValues' the
        // coefficients alpha can be used

        this->pointHierarchies[d]->getPoints(l, orderingConfiguration[d]);
        if (needsParam) {
          eval->setParameter(parameters[paramIndex]);
        }
        currentEvaluators.push_back(eval);
      }
      basisValues[d] = currentEvaluators[currentLevel]->getBasisValues();
      multiBounds[d] = this->pointHierarchies[d]->getNumPoints(currentLevel);

      if (needsParam) {
        ++paramIndex;
      }
    }

    updateBasisValues = false;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
