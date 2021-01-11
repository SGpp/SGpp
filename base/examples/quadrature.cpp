// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_quadrature_cpp Quadrature in C++
 *
 * The following example shows how to integrate in SG++, using both
 * direct integration of a sparse grid function and the use of
 * Monte Carlo integration.
 *
 * As in the \ref example_tutorial_cpp example, we deal with the function
 * \f[
 *   f\colon [0, 1]^2 \to \mathbb{R},\quad
 *   f(x_0, x_1) := 16 (x_0 - 1) x_0 (x_1 - 1) x_1
 * \f]
 * which we first interpolate. We then integrate the interpolant, then
 * the function itself using 100000 Monte Carlo points, and we then
 * compute the L2-error.
 *
 * The function, which sgpp::base::OperationQuadratureMC takes, has three parameters.
 * First, the dimensionality (int),
 * then a double* with the coordinates of the grid point \f$\in[0,1]^d\f$,
 * and finally a void* with clientdata for the function, see \ref sgpp::base::FUNC.
 */

// include all SG++ base headers
#include <sgpp_base.hpp>

#include <iostream>

using sgpp::base::OperationHierarchisation;

// function to interpolate
double f(int dim, double* x, void* clientdata) {
  double res = 1.0;

  for (int i = 0; i < dim; i++) {
    res *= 4.0 * x[i] * (1.0 - x[i]);
  }

  return res;
}

int main() {
  /**
   * Create a two-dimensional piecewise bi-linear grid of level 3
   */
  int dim = 2;
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  int level = 3;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;

  /**
   * Calculate the surplus vector alpha for the interpolant of \f$
   * f(x)\f$.  Since the function can be evaluated at any
   * point. Hence. we simply evaluate it at the coordinates of the
   * grid points to obtain the nodal values. Then we use
   * hierarchization to obtain the surplus value.
   *
   */

  sgpp::base::DataVector alpha(gridStorage.getSize());
  double p[2];

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    p[0] = gp.getStandardCoordinate(0);
    p[1] = gp.getStandardCoordinate(1);
    alpha[i] = f(2, p, nullptr);
  }

  std::unique_ptr<OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid))
      ->doHierarchisation(alpha);

  /**
     * Now we compute and compare the quadrature using four different methods available in SG++.
     */

  // direct quadrature
  std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(alpha);
  std::cout << "exact integral value:  " << res << std::endl;

  // Monte Carlo quadrature using 100000 paths
  sgpp::base::OperationQuadratureMC opMC(*grid, 100000);
  res = opMC.doQuadrature(alpha);
  std::cout << "Monte Carlo value:     " << res << std::endl;
  res = opMC.doQuadrature(alpha);
  std::cout << "Monte Carlo value:     " << res << std::endl;

  // Monte Carlo quadrature of a function
  res = opMC.doQuadratureFunc(f, nullptr);
  std::cout << "MC value:              " << res << std::endl;

  // Monte Carlo quadrature of error
  res = opMC.doQuadratureL2Error(f, nullptr, alpha);
  std::cout << "MC L2-error:           " << res << std::endl;
}  // end of main

/**
 * This results in an output similar to:
 * \verbinclude quadrature.output.txt
 */
