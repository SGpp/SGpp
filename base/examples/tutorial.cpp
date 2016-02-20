// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// All SG++ headers
// #include <sgpp_base.hpp>

// Or, better!, include only those that are required
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <iostream>

using SGPP::base::DataVector;
using SGPP::base::Grid;
using SGPP::base::GridGenerator;
using SGPP::base::GridIndex;
using SGPP::base::GridStorage;
using SGPP::base::OperationEval;

// function to interpolate
SGPP::float_t f(SGPP::float_t x0, SGPP::float_t x1) {
  return 16.0 * (x0 - 1) * x0 * (x1 - 1) * x1;
}

int main() {
  // create a two-dimensional piecewise bilinear grid
  size_t dim = 2;
  std::unique_ptr<Grid> grid = Grid::createLinearGrid(dim);
  GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:         " << gridStorage->dim() << std::endl;

  // create regular grid, level 3
  size_t level = 3;
  GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points:  " << gridStorage->size() << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage->size());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector: " << alpha.getSize() << std::endl;

  // set function values in alpha
  GridIndex* gp;

  for (size_t i = 0; i < gridStorage->size(); i++) {
    gp = gridStorage->get(i);
    alpha[i] = f(gp->getCoord(0), gp->getCoord(1));
  }

  std::cout << "alpha before hierarchization: " << alpha.toString() << std::endl;

  // hierarchize
  SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(
    alpha);
  std::cout << "alpha after hierarchization:  " << alpha.toString() << std::endl;

  // evaluate
  DataVector p(dim);
  p[0] = 0.52;
  p[1] = 0.73;
  OperationEval* opEval = SGPP::op_factory::createOperationEval(*grid);
  std::cout << "u(0.52, 0.73) = " << opEval->eval(alpha, p) << std::endl;
}
