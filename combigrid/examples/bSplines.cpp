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
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <sgpp/combigrid/operation/multidim/fullgrid/SLECoefficientsStorage.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <vector>

// Marsden: f(x) = 1 => coefficients(at least interior ones) are all 1
double f1D(sgpp::base::DataVector const &v) { return 1; }
double f2D(sgpp::base::DataVector const &v) { return v[0] * v[1]; }

int main() {
  size_t d = 1;
  sgpp::combigrid::MultiFunction func(f1D);

  sgpp::combigrid::CombiHierarchies::Collection grids(
      d, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      d, sgpp::combigrid::CombiEvaluators::BSplineInterpolation());
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());  // TODO(rehmemk): choose one

  // stores the values of the objective function
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
    // We store the results (= coefficients for Bspline interpolation) for each grid point, encoded
    // by a MultiIndex, in a TreeStorage
    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double> >(d);

    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    sgpp::combigrid::CombiEvaluators::Collection evalCopy(d);
    for (size_t dim = 0; dim < d; ++dim) {
      evalCopy[dim] = evaluators[dim]->cloneLinear();
      bool needsSorted = true;
      auto GridPoints = grids[dim]->getPoints(grid->getLevel()[dim], needsSorted);
      evalCopy[dim]->setGridPoints(GridPoints);
    }
    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields all multi-indices of grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
    auto it2 = funcStorage->getGuidedIterator(grid->getLevel(), it, std::vector<bool>(d, true));

    for (size_t index1 = 0; it2->isValid(); ++index1, it2->moveToNext()) {
      // Customize this computation for your algorithm
      double functionValue = it2->value();

      functionValues[index1] = functionValue;
      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      std::vector<std::vector<double> > basisValues;

      auto gridPoint = grid->getGridPoint(it.getMultiIndex());
      for (size_t dim = 0; dim < d; ++dim) {
        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
        auto basisValues1D = evalCopy[dim]->getBasisValues();
        std::vector<double> basisValues1D_vec(basisValues1D.size());
        for (size_t i = 0; i < basisValues1D.size(); i++) {
          basisValues1D_vec[i] = basisValues1D[i].value();
        }
        basisValues.push_back(basisValues1D_vec);  // basis values at grid points
      }

      for (size_t index2 = 0; innerIter.isValid(); ++index2, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < d; ++dim) {
          splineValue *= basisValues[dim][innerIndex[dim]];
        }
        A.set(index1, index2, splineValue);
      }
    }

    //    std::cout << A.toString() << std::endl;

    sgpp::optimization::FullSLE sle(A);
    sgpp::optimization::sle_solver::Armadillo solver;
    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    if (!solved) {
      exit(-1);
    }

    // ToDo (rehmemk) Linearkomb coefficients_sle und basisValues bilden

    it.reset();

    //    std::cout << A.size() << std::endl;
    //    for (size_t i = 0; i < coefficients_sle.size(); i++) {
    //      std::cout << coefficients_sle[i] << " ";
    //    }
    //    std::cout << "\n";

    // Ergebnis herausschreiben
    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
      coefficientTree->set(it.getMultiIndex(), coefficients_sle[vecIndex]);
    }

    return coefficientTree;
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
  parameter.set(0, 0.5);
  parameter.set(1, 0.5);
  parameter.set(2, 0.9);

  size_t levels = 3;
  double result = operation->evaluate(levels, parameter);

  std::cout << "Target function value: " << func(parameter) << "\n";
  std::cout << "Numerical result: " << result << "\n";
}
