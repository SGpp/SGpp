// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// All SG++ base headers
#include <sgpp_base.hpp>

#include <string>
#include <iostream>

using sgpp::base::Grid;
using sgpp::base::DataVector;
using sgpp::base::GridGenerator;
using sgpp::base::GridIndex;
using sgpp::base::GridStorage;
using sgpp::base::OperationEval;

int main() {
  //  size_t dim = 3;
  //  std::unique_ptr<Grid> grid = Grid::createPolyGrid(dim, 10);
  //  grid->getGenerator().regular(5);
  //
  //  std::filebuf fb;
  //  fb.open("new_poly.grid", std::ios::out);
  //  std::ostream os(&fb);
  //  os << grid->serialize();
  //  fb.close();

  std::ifstream ifs("poly.grid");
  std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
  std::unique_ptr<Grid> new_grid = Grid::unserialize(content);
  std::cout << new_grid->getSize() << std::endl;
}
