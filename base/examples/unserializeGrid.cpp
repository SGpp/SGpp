// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page unserializeGrid.cpp Grid unserialization Example
 *
 * In this example we show how to store a grid into a file and how to load it back
 * into a sgpp::base::Grid object.
 */

// We only need the header for dealing with the grid itself here
#include <sgpp/base/grid/Grid.hpp>

#include <string>
#include <iostream>

using sgpp::base::Grid;

int main() {
  /**
   * First, we create a four-dimensional linear grid of level 5
   */
  size_t dim = 4;
  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  grid->getGenerator().regular(5);
  std::cout << "Size of grid before serialization: " << grid->getSize() << std::endl;

  /**
    * Next, we store the grid to into the /tmp directory
    */
  std::filebuf fb;
  fb.open("/tmp/sgde-grid-4391dc6e-54cd-4ca2-9510-a9c02a2889ec.grid", std::ios::out);
  std::ostream os(&fb);
  os << grid->serialize();
  fb.close();

  /**
    * At last, we load the grid from the file back into a new object
    */
  std::ifstream ifs("/tmp/sgde-grid-4391dc6e-54cd-4ca2-9510-a9c02a2889ec.grid");
  std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
  std::unique_ptr<Grid> new_grid(Grid::unserialize(content));
  std::cout << "Size of grid after unserialization: " << new_grid->getSize() << std::endl;
}
