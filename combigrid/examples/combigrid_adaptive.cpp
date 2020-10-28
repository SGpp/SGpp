// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_combigrid_adaptive_cpp Combigrid Example Dimensional Adaptivity (C++)
 *
 * In this example, we use the combigrid module to adapt a combination grid solution to
 * best interpolate a test function at a given point.
 *
 * First, we include the required modules.
 */

#include <sgpp_base.hpp>
#include <sgpp_combigrid.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace std {
// needed for printing vectors
using sgpp::base::operator<<;
}  // namespace std

/**
 * We define parameters and perform hierarchization exactly as in the combigrid example
 */
int main() {
  // dimensionality
  const size_t dim = 2;
  // regular level
  const size_t n = 4;
  // B-spline degree
  const size_t p = 3;
  // whether there are points on the boundary
  const bool hasBoundary = true;
  // test function
  auto f = [](const sgpp::base::DataVector& x) {
    return std::sin(7.0 * x[0] - 3.0) * std::cos(5.0 * x[1] - 5.0);
  };

  // disable log output
  sgpp::base::Printer::getInstance().setVerbosity(-1);

  sgpp::base::SBsplineBase basis1d(p);
  sgpp::combigrid::HeterogeneousBasis basis(dim, basis1d);
  sgpp::combigrid::CombinationGrid combiGrid =
      sgpp::combigrid::CombinationGrid::fromRegularSparse(dim, n, basis, hasBoundary);

  sgpp::base::HashGridStorage gridStorage(dim);
  combiGrid.combinePoints(gridStorage);

  sgpp::base::DataVector fX(gridStorage.getSize());
  sgpp::base::DataVector x(dim);
  for (size_t k = 0; k < gridStorage.getSize(); k++) {
    gridStorage.getPoint(k).getStandardCoordinates(x);
    fX[k] = f(x);
  }

  std::vector<sgpp::base::DataVector> values;
  combiGrid.distributeValuesToFullGrids(gridStorage, fX, values);

  std::vector<sgpp::base::DataVector> surpluses(values);
  std::vector<std::unique_ptr<sgpp::combigrid::OperationPole>> opPole;
  sgpp::combigrid::OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(basis, opPole);
  sgpp::combigrid::OperationUPCombinationGrid opHier(combiGrid, opPole);
  opHier.apply(surpluses);

  /**
   * Let's assume that this time we are even more interested in the interpolated value
   * at a given point x, so much so that we would like to extend the combination scheme
   * in a way that makes this value more accurate.
   * f(x) can be considered the quantity of interest or QoI.
   */
  x.assign({0.12, 0.34});
  std::cout << "Value of test function at [" << x[0] << " " << x[1] << "]: " << f(x) << "\n";

  // create operation for evaluating and evaluate
  {
    sgpp::combigrid::OperationEvalCombinationGrid opEval(combiGrid);
    const double y = opEval.eval(surpluses, x);
    std::cout << "Value of combined sparse grid interpolant at [" << x[0] << " " << x[1]
              << "]: " << y << "\n";
  }

  /**
   * We start an AdaptiveCombinationGridGenerator from the full grids that are already in the
   * combiGrid we have now. This object is going to store the LevelVectors that are currently
   * adapted and the QoI results we know about them, as well as QoIs of grids that are not currently
   * adapted to yet.
   */

  // default parameters are linear summation for the QoIs, AveragingPriorityEstimator, and
  // WeightedRelevanceCalculator
  sgpp::combigrid::AdaptiveCombinationGridGenerator adaptiveCombinationGridGenerator =
      sgpp::combigrid::AdaptiveCombinationGridGenerator::fromCombinationGrid(combiGrid);

  // for each of the full grids, we calculate the value at point x and hand it to the generator
  for (size_t fullGridIndex = 0; fullGridIndex < combiGrid.getFullGrids().size(); ++fullGridIndex) {
    const sgpp::combigrid::FullGrid& fullGrid = combiGrid.getFullGrids()[fullGridIndex];
    const sgpp::combigrid::LevelVector& l = fullGrid.getLevel();
    {
      sgpp::combigrid::OperationEvalFullGrid opEval(fullGrid);
      const double y = opEval.eval(surpluses[fullGridIndex], x);
      std::cout << "Value of full grid [" << l[0] << " " << l[1] << "] interpolant at [" << x[0]
                << " " << x[1] << "]: " << y << "\n";
      adaptiveCombinationGridGenerator.setQoIInformation(l, y);
    }
  }

  // for (auto& l : adaptiveCombinationGridGenerator.getLevels()) {
  //   std::cout << l << ", ";
  // }
  // std::cout << std::endl;

  // all of the full grids known so far should already be in the generator's 'old set', so this call
  // here would return false:
  bool adapted = adaptiveCombinationGridGenerator.adaptAllKnown();
  assert(!adapted);

  // looking at the 'active set', we see which full grid spaces could potentially be interesting
  // next, they are the upper neighbors of the old set
  const auto activeSet = adaptiveCombinationGridGenerator.getActiveSet();

  // we also get the values at x for those full grids and add the QoI information to the generator
  for (const auto& levelVector : activeSet) {
    const sgpp::combigrid::FullGrid fullGrid{levelVector, basis, hasBoundary};

    // to evaluate, we fill the full grid with values of f
    // TODO(pollinta) this feels like a really ugly way to evaluate a function on the grid
    // TODO(pollinta) I am happy about syntax suggestions to make this more readable
    sgpp::base::DataVector fullGridValues(fullGrid.getNumberOfIndexVectors());
    auto range = sgpp::combigrid::IndexVectorRange(fullGrid);
    for (auto it = range.begin(); it != range.end(); ++it) {
      sgpp::base::DataVector z(dim);
      it.getStandardCoordinates(z);
      fullGridValues[it.getSequenceNumber()] = f(z);
      // std::cout << "Value of full grid [" << levelVector[0] << " " << levelVector[1] << "] at ["
      //           << z[0] << " " << z[1] << "]: " << f(z) << "\n";
    }

    // like before, hierarchize the full grid to evaluate using the basis functions
    sgpp::combigrid::OperationUPCombinationGrid opHierSingleGrid(
        sgpp::combigrid::CombinationGrid(fullGrid), opPole);
    auto fullGridSurpluses = std::vector<sgpp::base::DataVector>(1, fullGridValues);
    opHierSingleGrid.apply(fullGridSurpluses);

    sgpp::combigrid::OperationEvalFullGrid opEval(fullGrid);
    const double y = opEval.eval(fullGridSurpluses[0], x);
    std::cout << "Value of new full grid [" << levelVector[0] << " " << levelVector[1]
              << "] interpolant at [" << x[0] << " " << x[1] << "]: " << y << "\n";
    adaptiveCombinationGridGenerator.setQoIInformation(levelVector, y);
  }

  for (auto& l : adaptiveCombinationGridGenerator.getOldSet()) {
    std::cout << l << ", ";
  }
  std::cout << std::endl;

  // now, we can adapt some grids
  adaptiveCombinationGridGenerator.adaptNextLevelVector();
  adaptiveCombinationGridGenerator.adaptNextLevelVector();

  for (auto& l : adaptiveCombinationGridGenerator.getOldSet()) {
    std::cout << l << ", ";
  }
  std::cout << std::endl;

  // if the computation would become expensive, e.g. for higher resolutions, we could get a priority
  // estimate which levels should be calculated first
  auto priorityQueue = adaptiveCombinationGridGenerator.getPriorityQueue();

  // finally, we may adapt the generator to all calculated values
  adaptiveCombinationGridGenerator.adaptAllKnown();

  /**
   * Now, we can generate a new combination grid that contains all the subspaces / full grids that
   * we have adapted to. We again interpolate the value at x and see whether it has become more
   * accurate.
   */
  auto combiGridAdapted = adaptiveCombinationGridGenerator.getCombinationGrid(basis, hasBoundary);

  // fill the combination grid with analytical values
  sgpp::base::HashGridStorage gridStorageAdapted(dim);
  combiGridAdapted.combinePoints(gridStorageAdapted);

  sgpp::base::DataVector fAdapted(gridStorageAdapted.getSize());
  for (size_t k = 0; k < gridStorageAdapted.getSize(); k++) {
    sgpp::base::DataVector z(dim);
    gridStorageAdapted.getPoint(k).getStandardCoordinates(z);
    fAdapted[k] = f(z);
  }

  // distribute the values to the full grids again and hierarchize
  std::vector<sgpp::base::DataVector> valuesAdapted;
  combiGridAdapted.distributeValuesToFullGrids(gridStorageAdapted, fAdapted, valuesAdapted);

  std::vector<sgpp::base::DataVector> surplusesAdapted(valuesAdapted);
  sgpp::combigrid::OperationUPCombinationGrid opHierAdapted(combiGridAdapted, opPole);
  opHierAdapted.apply(surplusesAdapted);

  {
    sgpp::combigrid::OperationEvalCombinationGrid opEval(combiGridAdapted);
    const double y = opEval.eval(surplusesAdapted, x);
    std::cout << "Value of combined sparse grid interpolant at [" << x[0] << " " << x[1]
              << "]: " << y << "\n";
  }

  return 0;
}

/**
 * The example program outputs the following results:
 * \verbinclude combigrid_adaptive.output.txt
 *
 * We see that the value of the adapted combined sparse grid interpolant at the
 * evaluation point is even
 * closer to the actual value of the test function than the smaller non-adapted combination grid we
 * used in the last example.
 */
