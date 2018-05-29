// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <iostream>

/**
 * The objective function we want to interpolate
 */
double objectiveFunction(sgpp::base::DataVector v) {
  double res = 1;
  for (size_t d = 0; d < v.getSize(); d++) {
    res *= v[d];
  }
  return res;
}

int main() {
  /**
   * First, we create a two-dimensional grid (type sgpp::base::Grid)
   * with B-spline basis functions
   */
  size_t dim = 1;
  size_t degree = 3;
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createBsplineBoundaryGrid(dim, degree));

  /**
   * Then we obtain a reference to the grid's
   * sgpp::base::GridStorage object which allows us, e.g., to access grid
   * points, to obtain the dimensionality (which we print) and the
   * number of grid points.
   */
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:         " << gridStorage.getDimension() << std::endl;

  /**
   * Now, we use a sgpp::base::GridGenerator to create a regular sparse grid of level 3 and print
   * its number of gridpoints.
   */
  size_t level = 4;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points:  " << gridStorage.getSize() << std::endl;

  /**
   * We create an object of type sgpp::base::DataVector
   * which is essentially a wrapper around a double array.
   * The DataVector is initialized with as many
   * entries as there are grid points. It serves as a coefficient vector for the
   * sparse grid interpolant we want to construct. As the entries of a
   * freshly created DataVector are not initialized, we set them to
   * 0.0. (This is superfluous here as we initialize them in the
   * next few lines anyway.)
   */
  sgpp::base::DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector: " << alpha.getSize() << std::endl;

  /**
   * evaluate the ibjective function in the grid points to obtain the right hand side of the system
   * of linear equations (SLE)
   */
  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    f_values[i] = objectiveFunction(p);
  }

  /**
   * Choose a solver type, initialize and solve the SLE
   */
  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }

  /**
   * Create the interpolant from the grid and the solution alpha of the above SLE
   */
  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*grid, alpha);

  /**
   * Evaluate the interpolant in (0.3,...0.3)
   */
  sgpp::base::DataVector evalPoint(dim, 0.3);
  std::cout << " f(p) = " << objectiveFunction(evalPoint)
            << " , I(p) = " << sparseGridSurrogate.eval(evalPoint) << " error = "
            << fabs(objectiveFunction(evalPoint) - sparseGridSurrogate.eval(evalPoint))
            << std::endl;

  return 0;
}
