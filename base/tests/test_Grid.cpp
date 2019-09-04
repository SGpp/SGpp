// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>

using sgpp::base::Grid;

BOOST_AUTO_TEST_SUITE(test_Grid)

BOOST_AUTO_TEST_CASE(test_clone) {
  size_t dim = 2;
  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  size_t level = 3;
  grid->getGenerator().regular(level);

  auto newGrid = grid->clone();

  BOOST_CHECK_EQUAL(newGrid->getStorage().getSize(), grid->getStorage().getSize());
  BOOST_CHECK_EQUAL(newGrid->getDimension(), grid->getDimension());

  delete newGrid;
}

BOOST_AUTO_TEST_SUITE_END()
