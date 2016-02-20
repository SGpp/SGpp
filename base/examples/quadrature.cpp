// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// All SG++ base headers
#include <sgpp_base.hpp>

#include <iostream>

using SGPP::base::DataVector;
using SGPP::base::Grid;
using SGPP::base::GridGenerator;
using SGPP::base::GridIndex;
using SGPP::base::GridStorage;
using SGPP::base::OperationQuadrature;
using SGPP::base::OperationQuadratureMC;

// function to interpolate
SGPP::float_t f(int dim, SGPP::float_t* x, void* clientdata) {
  SGPP::float_t res = 1.0;

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
  std::cout << "dimensionality:        " << gridStorage.dim() << std::endl;

  // create regular grid, level 3
  int level = 3;
  std::unique_ptr<GridGenerator> gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage.size() << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage.size());
  GridIndex* gp;
  SGPP::float_t p[2];

  for (size_t i = 0; i < gridStorage.size(); i++) {
    gp = gridStorage.get(i);
    p[0] = gp->getCoord(0);
    p[1] = gp->getCoord(1);
    alpha[i] = f(2, p, NULL);
  }

  SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(
    alpha);

  // direct quadrature
  OperationQuadrature* opQ = SGPP::op_factory::createOperationQuadrature(*grid);
  SGPP::float_t res = opQ->doQuadrature(alpha);
  std::cout << "exact integral value:  " << res << std::endl;
  delete opQ;

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
