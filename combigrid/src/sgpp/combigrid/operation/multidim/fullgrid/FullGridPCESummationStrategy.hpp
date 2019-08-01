// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This is an experimental class using an alternative basis change for PCE. Please use other classes
 * for standard PCE instead.
 */
template <typename V>
class FullGridPCESummationStrategy : public AbstractFullGridSummationStrategy<V> {
 public:
  FullGridPCESummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> evaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, evaluatorPrototypes, pointHierarchies) {
    size_t numDimensions = pointHierarchies.size();
    for (size_t d = 0; d < numDimensions; d++) {
      auto evalConfig = this->evaluatorPrototypes[d]->getConfig();
      auto basisType = evalConfig.type;
      if (basisType == CombiEvaluatorTypes::Tensor_PolynomialInterpolation) {
        functionBases.push_back(evalConfig.functionBasis);
      } else {
        throw sgpp::base::algorithm_exception(
            "FullGridPCESummationStrategy: this evaluator is currently "
            "not supported.");
      }
    }
  }

  virtual ~FullGridPCESummationStrategy() {}

  V eval(MultiIndex const &level) override {
    size_t numDimensions = this->pointHierarchies.size();
    std::vector<bool> orderingConfiguration(numDimensions, false);

    auto grid = this->getTensorGrid(level, orderingConfiguration);
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
    auto funcIter = this->storage->getGuidedIterator(level, it, orderingConfiguration);

    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
      auto gridPoint = grid->getGridPoint(funcIter->getMultiIndex());
      functionValues[ixEvalPoints] = funcIter->value();

      std::vector<std::vector<double>> basisValues;
      for (size_t dim = 0; dim < numDimensions; ++dim) {
        std::vector<double> basisValues1D_vec(numGridPointsVec[dim]);
        for (size_t i = 0; i < numGridPointsVec[dim]; i++) {
          basisValues1D_vec[i] = functionBases[dim]->evaluate(i, gridPoint[dim]);
        }
        basisValues.push_back(basisValues1D_vec);
      }

      sgpp::combigrid::MultiIndexIterator innerIter(numGridPointsVec);
      for (size_t ixBasisFunctions = 0; innerIter.isValid();
           ++ixBasisFunctions, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < numDimensions; ++dim) {
          splineValue *= basisValues[dim][innerIndex[dim]];
        }
        A.set(ixEvalPoints, ixBasisFunctions, splineValue);
      }
    }

    sgpp::base::FullSLE sle(A);
    sgpp::base::sle_solver::Auto solver;
    sgpp::base::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    if (!solved) {
      std::cout << "FullGridPCEEvaluator::eval(): system could not be solved\n";
    }

    it.reset();

    auto resultCoefficients =
        std::make_shared<sgpp::combigrid::TreeStorage<FloatScalarVector>>(numDimensions);
    for (size_t vecIndex_i = 0; it.isValid(); ++vecIndex_i, it.moveToNext()) {
      resultCoefficients->set(it.getMultiIndex(), FloatScalarVector(coefficients_sle[vecIndex_i]));
    }

    return V(resultCoefficients);
  }

 private:
  std::vector<std::shared_ptr<AbstractInfiniteFunctionBasis1D>> functionBases;

  /**
   * @return Returns the grid in the current level as a pointer to a TensorGrid object
   */
  std::shared_ptr<TensorGrid> getTensorGrid(MultiIndex const &level,
                                            std::vector<bool> orderingConfiguration) {
    size_t numDimensions = this->pointHierarchies.size();
    std::vector<base::DataVector> grids1D;

    for (size_t d = 0; d < numDimensions; ++d) {
      bool sorted = orderingConfiguration[d];
      grids1D.push_back(base::DataVector(this->pointHierarchies[d]->getPoints(level[d], sorted)));
    }

    return std::make_shared<TensorGrid>(grids1D, level);
  }
};

} /* namespace combigrid */
} /* namespace sgpp */
