// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


// Include all SGpp Base headers
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
    alpha[i] = f(2, p, NULL);
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
  res = opMC.doQuadratureFunc(f, NULL);
  std::cout << "MC value:              " << res << std::endl;

  // Monte Carlo quadrature of error
  res = opMC.doQuadratureL2Error(f, NULL, alpha);
  std::cout << "MC L2-error:           " << res << std::endl;
}
