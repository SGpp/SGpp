// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include "../../../../../../../base/src/sgpp/base/tools/Printer.hpp"
#include "../../../../../../../base/src/sgpp/base/tools/sle/solver/Auto.hpp"
#include "../../../../../../../base/src/sgpp/base/tools/sle/system/FullSLE.hpp"

#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#endif

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This is an experimental class using an alternative basis change for PCE. Please use other classes
 * for standard PCE instead. // TODO(holzmudd)
 */
template <typename V>
class FullGridOptimizedPCESummationStrategy : public AbstractFullGridSummationStrategy<V> {
 public:
  FullGridOptimizedPCESummationStrategy(
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<V>>> scalarProductEvaluatorPrototypes,
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies)
      : AbstractFullGridSummationStrategy<V>(storage, scalarProductEvaluatorPrototypes,
                                             pointHierarchies) {
#ifdef USE_EIGEN
    inverseMatrices.resize(pointHierarchies.size());
#endif
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

  virtual ~FullGridOptimizedPCESummationStrategy() {}

  V eval(MultiIndex const &level) override {
#ifdef USE_EIGEN
    size_t numDimensions = this->pointHierarchies.size();
    std::vector<bool> orderingConfiguration(numDimensions, false);

    auto grid = this->getTensorGrid(level, orderingConfiguration);
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }
    auto grids1D = grid->get1DGrids();  // TODO(holzmudd): evtl. optimierbar

    // compute 1D basis change matrices if not already computed
    for (size_t dim = 0; dim < numDimensions; ++dim) {
      while (inverseMatrices[dim].size() <= level[dim]) {
        size_t level1D = inverseMatrices[dim].size();
        size_t numPoints = this->pointHierarchies[dim]->getNumPoints(level1D);
        Eigen::MatrixXd mat(numPoints, numPoints);
        for (size_t i = 0; i < numPoints; ++i) {
          for (size_t j = 0; j < numPoints; ++j) {
            mat(i, j) =
                functionBases[dim]->evaluate(j, this->pointHierarchies[dim]->getPoint(level1D, i));
          }
        }
        inverseMatrices[dim].push_back(mat.fullPivHouseholderQr().inverse());
      }
    }

    auto resultCoefficients =
        std::make_shared<sgpp::combigrid::TreeStorage<FloatScalarVector>>(numDimensions);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
    auto funcIter = this->storage->getGuidedIterator(level, it, orderingConfiguration);

    std::vector<double> functionValues(numGridPoints);

    // store all function values into the vector for faster access (needed often)
    for (size_t i = 0; funcIter->isValid(); ++i, funcIter->moveToNext()) {
      functionValues[i] = funcIter->value();
    }

    it.reset();

    // for all relevant basis functions
    for (; it.isValid(); it.moveToNext()) {
      double sum = 0.0;  // sum over contributions from all function values

      sgpp::combigrid::MultiIndexIterator innerIt(grid->numPoints());
      // for all grid points
      for (size_t i = 0; innerIt.isValid(); ++i, innerIt.moveToNext()) {
        double prod = functionValues[i];
        for (size_t dim = 0; dim < numDimensions; ++dim) {
          prod *= inverseMatrices[dim][level[dim]](it.indexAt(dim), innerIt.indexAt(dim));
        }
        sum += prod;
      }

      resultCoefficients->set(it.getMultiIndex(), FloatScalarVector(sum));
    }

    return V(resultCoefficients);
#else
    throw sgpp::base::application_exception("need Eigen to use the PCE transformation.");
#endif
  }

 private:
  std::vector<std::shared_ptr<AbstractInfiniteFunctionBasis1D>> functionBases;
#ifdef USE_EIGEN
  std::vector<std::vector<Eigen::MatrixXd>> inverseMatrices;  // indices: dim, level
#endif

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
