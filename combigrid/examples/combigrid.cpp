// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_combigrid_cpp Combigrid Example (C++)
 *
 * In this example, we use the combigrid module to interpolate a test function on a two-dimensional
 * regular sparse grid with the combination technique and hierarchical B-splines.
 *
 * First, we include the required modules.
 */
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "sgpp_base.hpp"
#include "sgpp_combigrid.hpp"

/**
 * We define some parameters such as dimensionality and level of the regular sparse grid.
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
        return std::sin(7.0*x[0]-3.0) * std::cos(5.0*x[1]-5.0);
      };

  // disable log output
  sgpp::base::Printer::getInstance().setVerbosity(-1);

  /**
   * The basis functions are defined via an sgpp::combigrid::HeterogeneousBasis object. In contrast
   * to sgpp::base::Basis, this allows for different types of basis functions for the different
   * dimensions. However, for this example, we do not need this flexibility, so we use the same
   * basis function types for both dimensions.
   */
  sgpp::base::SBsplineBase basis1d(p);
  sgpp::combigrid::HeterogeneousBasis basis(dim, basis1d);

  /**
   * An sgpp::combigrid::CombinationGrid is a collection of full grids (nodal subspaces) together
   * with scalar-valued coefficients. Here, we construct an sgpp::combigrid::CombinationGrid object
   * for a regular sparse grid via the combination technique.
   */
  sgpp::combigrid::CombinationGrid combiGrid =
      sgpp::combigrid::CombinationGrid::fromRegularSparse(dim, n, basis, hasBoundary);

  /**
   * We obtain the grid points of the regular sparse grid by combining the grid points of all
   * full grids that are contained in the combination grid.
   */
  sgpp::base::HashGridStorage gridStorage(dim);
  combiGrid.combinePoints(gridStorage);

  // evaluate test function at grid points
  sgpp::base::DataVector fX(gridStorage.getSize());
  sgpp::base::DataVector x(dim);

  for (size_t k = 0; k < gridStorage.getSize(); k++) {
    gridStorage.getPoint(k).getStandardCoordinates(x);
    fX[k] = f(x);
  }

  /**
   * We now want to perform an operation on each full grid. For this, we distribute the values of
   * the combined grid (sparse grid) to the full grids. The result is a \c std::vector of
   * sgpp::base::DataVector; each \c DataVector contains the values at all grid points for one
   * specific full grid.
   */
  std::vector<sgpp::base::DataVector> values;
  combiGrid.distributeValuesToFullGrids(gridStorage, fX, values);

  /**
   * The operation we want to perform on each full grid is hierarchization. Since the grids are
   * full grids, we can use the unidirectional principle for this, which performs 1D hierarchization
   * on each pole (one-dimensional sub-grid), iterating over all dimensions.
   */
  // copy the values (surpluses will be modified in-place)
  std::vector<sgpp::base::DataVector> surpluses(values);

  // create pole operation
  std::vector<std::unique_ptr<sgpp::combigrid::OperationPole>> opPole;
  sgpp::combigrid::OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(basis, opPole);

  // create operation for unidirectional principle and hierarchize in-place
  sgpp::combigrid::OperationUPCombinationGrid opHier(combiGrid, opPole);
  opHier.apply(surpluses);

  /**
   * The resulting surpluses are also a \c std::vector of sgpp::base::DataVector, separated by full
   * grids. We could combine the full grid surpluses via the combination formula to the sparse grid
   * surpluses via \c combineSparseGridValues. However, the operation
   * sgpp::combigrid::OperationEvalCombinationGrid does this automatically.
   *
   * We evaluate the combined function (combination of all full grid interpolants) at some arbitrary
   * point and print the value.
   */
  // test point at which to evaluate
  x.assign({0.12, 0.34});
  std::cout << "Value of test function at [" << x[0] << " " << x[1] << "]: " << f(x) << "\n";

  // create operation for evaluating and evaluate
  {
    sgpp::combigrid::OperationEvalCombinationGrid opEval(combiGrid);
    const double y = opEval.eval(surpluses, x);
    std::cout << "Value of combined sparse grid interpolant at [" <<
        x[0] << " " << x[1] << "]: " << y << "\n";
  }

  /**
   * Finally, we do the same for one full grid of the combination grid: We evaluate the
   * corresponding interpolant. We extract the surpluses from the already calculated \c vector of
   * \c DataVector. Alternatively, we could also apply sgpp::combigrid::OperationUPFullGrid with
   * opPole to obtain the surpluses for this single full grid.
   */
  // select the second full grid of the combination grid (arbitrary choice)
  const size_t fullGridIndex = 1;
  const sgpp::combigrid::FullGrid& fullGrid = combiGrid.getFullGrids()[fullGridIndex];
  const sgpp::combigrid::LevelVector& l = fullGrid.getLevel();
  std::cout << "Level of selected full grid with index " << fullGridIndex << ": [" <<
      l[0] << " " << l[1] << "]\n";

  // create operation for evaluating and evaluate
  {
    sgpp::combigrid::OperationEvalFullGrid opEval(fullGrid);
    const double y = opEval.eval(surpluses[fullGridIndex], x);
    std::cout << "Value of full grid interpolant at [" << x[0] << " " << x[1] << "]: " << y << "\n";
  }

  return 0;
}

/**
 * The example program outputs the following results:
 * \verbinclude combigrid.output.txt
 *
 * We see that the value of the combined sparse grid interpolant at the evaluation point is
 * closer to the actual value of the test function than the value of the chosen full grid
 * interpolant, which corresponds to the full grid of level \f$(3, 1)\f$.
 */
