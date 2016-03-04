// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// All SG++ base headers
#include <sgpp_base.hpp>

#include <iostream>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridIndex;
using sgpp::base::GridStorage;
using sgpp::base::OperationQuadrature;
using sgpp::base::OperationQuadratureMC;

// function to interpolate
double f(int dim, double* x, void* clientdata) {
  double res = 1.0;

  for (int i = 0; i < dim; i++) {
    res *= 4.0 * x[i] * (1.0 - x[i]);
  }

  return res;
}

int main() {
  // create a two-dimensional piecewise bi-linear grid
  int dim = 2;
  std::unique_ptr<Grid> grid = Grid::createLinearGrid(dim);
  GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  int level = 3;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());
  GridIndex* gp;
  double p[2];

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    gp = gridStorage.get(i);
    p[0] = gp->getCoord(0);
    p[1] = gp->getCoord(1);
    alpha[i] = f(2, p, NULL);
  }

  sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(
    alpha);

  // direct quadrature
  std::unique_ptr<OperationQuadrature> opQ(sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(alpha);
  std::cout << "exact integral value:  " << res << std::endl;

  // Monte Carlo quadrature using 100000 paths
  OperationQuadratureMC opMC(*grid, 100000);
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
