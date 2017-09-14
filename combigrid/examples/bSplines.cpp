// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const &v) { return v[0]; }

/**
 * @section combigrid_example_6 Example 6: Using a function operating on grids
 *
 * This example shows how to apply different operators in different dimensions.
 */
int main() {
  size_t d = 3;

  sgpp::combigrid::MultiFunction func(f);

  /**
   * To create a CombigridOperation, we currently have to use the longer way as in example 5.
   */
  sgpp::combigrid::CombiHierarchies::Collection grids(
      d, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      d, sgpp::combigrid::CombiEvaluators::
             cubicSplineInterpolation());  // TODO(rehmemk): bSplineInterpolation(degree)
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());  // TODO(rehmemk): choose one

  auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

  // sgpp::combigrid::MultiFunction dummyFunc([](sgpp::base::DataVector const &v) { return 0.0; });

  // auto bSplineCoefficientsStorage =
  //     std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, false, dummyFunc);

  /**
   * In some applications, you might not want to have a callback function that is called at single
   * points, but on a full grid. One of these applications is solving PDEs. This example provides a
   * simple framework where a PDE solver can be included. It is also suited for other tasks.
   * The core part is a function that computes grid values on a full grid.
   */
  sgpp::combigrid::GridFunction gf([&](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    // We store the results for each grid point, encoded by a MultiIndex, in a TreeStorage
    auto result = std::make_shared<sgpp::combigrid::TreeStorage<double> >(d);

    sgpp::combigrid::CombiEvaluators::Collection evalCopy(d);
    for (size_t dim = 0; dim < d; ++dim) {
      evalCopy[dim] = evaluators[dim]->cloneLinear();
      bool needsSorted = true;
      evalCopy[dim]->setGridPoints(grids[dim]->getPoints(grid->getLevel()[dim], needsSorted));
    }

    // Creates an iterator that yields all multi-indices of grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
    auto it2 = funcStorage->getGuidedIterator(grid->getLevel(), it, std::vector<bool>(d, true));

    for (size_t index1 = 0; it2->isValid(); ++index1, it2->moveToNext()) {
      // Customize this computation for your algorithm
      double functionValue = it2->value();

      // functionValues(index1) = functionValue

      // Matrix aufstellen

      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      std::vector<std::vector<double> > basisCoefficients;

      auto gridPoint = grid->getGridPoint(it.getMultiIndex());
      for (size_t dim = 0; dim < d; ++dim) {
        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
        auto basisValues = evalCopy[dim]->getBasisValues();
        std::vector<double> basisValues_vec(basisValues.size());
        for (size_t i = 0; i < basisValues.size(); i++) {
          basisValues_vec[i] = basisValues[i].value();
        }
        basisCoefficients.push_back(basisValues_vec);  // basis values at grid points
      }

      for (size_t index2 = 0; innerIter.isValid(); ++index2, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < d; ++dim) {
          splineValue *= basisCoefficients[dim][innerIndex[dim]];
        }

        // matrix(index1, index2) = splineValue
      }
    }

    // solution = matrix.inv() * functionValues

    // Gleichungssystem lösen

    it.reset();

    // Ergebnis herausschreiben
    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
      result->set(it.getMultiIndex(), 5.0);  // solution(vecIndex)
    }

    return result;
  });

  /**
   * We have to specify if the function always produces the same value for the same grid points.
   * This can make the storage smaller if the grid points are nested. In this implementation, this
   * is true. However, it would be false in the PDE case, so we set it to false here.
   */
  bool exploitNesting = false;

  /**
   * Now create an operation as usual and evaluate the interpolation with a test parameter.
   */
  auto operation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      grids, evaluators, levelManager, gf, exploitNesting);

  sgpp::base::DataVector parameter(d);
  parameter.set(0, 0.1);
  parameter.set(1, 0.2);
  parameter.set(2, 0.3);

  double result = operation->evaluate(4, parameter);

  std::cout << "Target function value: " << func(parameter) << "\n";
  std::cout << "Numerical result: " << result << "\n";
}
