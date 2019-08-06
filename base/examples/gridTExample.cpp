// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_gridTExample_cpp Generalised Sparse Grids 
 * This example creates a generalised grid. The first CLI argument is
 * the number of dimensions, the second the level and the third the parameter
 * T. It then prints out the grid size.
 */

#include <sgpp/base/grid/Grid.hpp>
#include <cstdlib>
#include <iostream>

int main(int argc, char **argv) {
  auto dimensions = ((argc >= 2) ? std::atoi(argv[1]) : 2);
  auto level = ((argc >= 3) ? std::atoi(argv[2]) : 3);
  double T = ((argc >= 4) ? std::atof(argv[3]) : 0);

  auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
  auto& generator = grid->getGenerator();
  generator.regular(level, T);
  std::cout << grid->getSize() << std::endl;
}
