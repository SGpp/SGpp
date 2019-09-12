// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_SCALAPACK
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <mpi.h>
#include <cmath>
#include <memory>

using sgpp::datadriven::BlacsProcessGrid;

BOOST_AUTO_TEST_SUITE(testBlacsProcessGrid)

BOOST_AUTO_TEST_CASE(testConstructor) {
  std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(0, 0);

  int availableProcesses = BlacsProcessGrid::availableProcesses();

  if (grid->isProcessInGrid()) {
    BOOST_CHECK_EQUAL(grid->getTotalRows(), static_cast<int>(std::sqrt(availableProcesses)));
    BOOST_CHECK_EQUAL(grid->getTotalColumns(), static_cast<int>(std::sqrt(availableProcesses)));
  }

  if (availableProcesses >= 2) {
    grid = std::make_shared<BlacsProcessGrid>(2, 1);

    if (grid->isProcessInGrid()) {
      BOOST_CHECK_EQUAL(grid->getTotalRows(), 2);
      BOOST_CHECK_EQUAL(grid->getTotalColumns(), 1);
    }

    grid = std::make_shared<BlacsProcessGrid>(1, 1);

    if (BlacsProcessGrid::getCurrentProcess() == 0) {
      BOOST_CHECK(grid->isProcessInGrid());
    } else {
      BOOST_CHECK(!grid->isProcessInGrid());
    }
  }
}

BOOST_AUTO_TEST_CASE(testNumberOfProcesses) {
  // test availableProcesses()
  int mpiTasks = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpiTasks);

  int availableProcesses = BlacsProcessGrid::availableProcesses();
  BOOST_CHECK_EQUAL(availableProcesses, mpiTasks);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SCALAPACK */
