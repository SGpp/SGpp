// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
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
class FullGridQuadraticSummationStrategy : public AbstractFullGridSummationStrategy<V> {
 protected:
 private:
  /*
   * vector product v * transpose(v2)
   */

  double dotMult(std::vector<double> const &vector1, std::vector<double> const &vector2) {
    double result(0.0);
    if (vector1.size() == vector2.size()) {
      for (size_t i = 0; i < vector1.size(); i++) {
        result += vector1[i] * vector2[i];
      }
    } else {
      std::cerr << "dotMult for vectors of different size is undefined. returning 0" << std::endl;
    }
    return result;
  }

  FloatScalarVector dotMult(FloatArrayVector vector1, FloatArrayVector vector2) {
    FloatScalarVector result(0.0);
    if (vector1.size() == vector2.size()) {
      auto values1 = vector1.getValues();
      auto values2 = vector2.getValues();
      size_t index = 0;
      for (auto &val : values1) {
        val.componentwiseMult(values2[index]);
        result.add(val);
        index++;
      }
    } else {
      //      std::cout << vector1.size() << " " << vector2.size() << std::endl;
      //      for (size_t i = 0; i < vector1.size(); i++) {
      //        std::cout << vector1[i] << " ";
      //      }
      //      std::cout << "\n";
      //      for (size_t i = 0; i < vector2.size(); i++) {
      //        std::cout << vector2[i] << " ";
      //      }
      //      std::cout << "\n";

      std::cerr << "dotMult for FloatArrayVectors of different size is undefined. returning 0"
                << std::endl;
    }
    return result;
  }

  FloatScalarVector dotMult(FloatScalarVector vector1, FloatArrayVector vector2) {
    std::cerr << "FullGridQuadraticSummationstrategy: operation dotMut is undefined for "
                 "FloatScalarVector."
              << std::endl;
    return FloatScalarVector::zero();
  }

  FloatScalarVector dotMult(FloatTensorVector vector1, FloatArrayVector vector2) {
    std::cerr << "FullGridQuadraticSummationstrategy: operation dotMut is undefined for "
                 "FloatScalarVector."
              << std::endl;
    return FloatScalarVector::zero();
  }

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

  ~FullGridQuadraticSummationStrategy() {}

  /**
   * Evaluates the function given through the storage for a certain level-multi-index (see class
   * description).
   * Summation of the form \sum_i \alpha_i \sum_j \alpha_j  basis_i(param) basis_j(param)
   *                    <=> v^T * A * v, A_ij = < basis_i(param),basis_j(param) >, v_i = alpha_i
   * This is used for variance calculations.
   * Currently it can only be used with template type V = FloatArrayVector
   */
  V eval(MultiIndex const &level) override {
    //    std::cout << "level " << level[0] << " " << level[1] << " " << level[2] << std::endl;

    CGLOG("FullGridTensorEvaluator::eval(): start");
    size_t numDimensions = this->evaluators.size();
    //    size_t lastDim = numDimensions - 1;
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

    CGLOG("FullGridTensorEvaluator::eval(): create storage iterator");

    // start iteration
    MultiIndexIterator it(multiBounds);

    std::shared_ptr<AbstractMultiStorageIterator<double>> funcIter =
        this->storage->getGuidedIterator(level, it, orderingConfiguration);
    CGLOG("FullGridTensorEvaluator::eval(): start loop");

    if (!funcIter->isValid()) {  // should not happen
      return V::zero();
    }

    std::vector<double> coefficients;
    while (true) {
      coefficients.push_back(funcIter->value());
      int h = funcIter->moveToNext();
      if (h < 0) {
        break;
      }
    }
    it.reset();
    std::vector<std::vector<double>> M;
    while (it.isValid()) {
      MultiIndexIterator innerIt(multiBounds);
      MultiIndex indexI = it.getMultiIndex();
      std::vector<double> lineI;
      while (innerIt.isValid()) {
        double entry = 1;
        MultiIndex indexJ = innerIt.getMultiIndex();
        for (size_t d = 0; d < numDimensions; d++) {
          entry *= this->basisValues[d][indexI[d]][indexJ[d]].value();
        }
        lineI.push_back(entry);
        innerIt.moveToNext();
      }
      M.push_back(lineI);
      it.moveToNext();
    }

    std::vector<double> innerSum;
    for (size_t i = 0; i < coefficients.size(); i++) {
      //      std::cout << M[i].size() << " " << coefficients.size() << std::endl;
      innerSum.push_back(dotMult(M[i], coefficients));
    }

    double variance = dotMult(innerSum, coefficients);
    FloatScalarVector var(variance);
    V result(var);
    return result;
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
