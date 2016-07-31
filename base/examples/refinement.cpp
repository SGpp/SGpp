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
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <iostream>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::SurplusRefinementFunctor;

// function to interpolate
double f(double x0, double x1) {
  return 16.0 * (x0 - 1) * x0 * (x1 - 1) * x1;
}

int main() {
  // create a two-dimensional piecewise bilinear grid
  size_t dim = 2;
  std::unique_ptr<Grid> grid = Grid::createLinearGrid(dim);
  GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:                   " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  size_t level = 3;
  grid->getGenerator().regular(level);
  std::cout << "number of initial grid points:    " << gridStorage.getSize() << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector:           " << alpha.getSize() << std::endl;

  // refine adaptively 5 times
  for (int step = 0; step < 5; step++) {
    // set function values in alpha
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      GridPoint& gp = gridStorage.getPoint(i);
      alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1));
    }

    // hierarchize
    sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(
      alpha);

    // refine a single grid point each time
    SurplusRefinementFunctor functor(alpha, 1);
    grid->getGenerator().refine(functor);
    std::cout << "refinement step " << step + 1 << ", new grid size: " << alpha.getSize()
         << std::endl;

    // extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.getSize());
  }
}
