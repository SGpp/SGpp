// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>

#include <boost/test/unit_test.hpp>

#include <climits>
#include <cmath>
#include <vector>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::HashCoarsening;
using sgpp::base::HashGenerator;
using sgpp::base::HashGridPoint;
using sgpp::base::HashGridStorage;
using sgpp::base::SurplusCoarseningFunctor;

/*
 * dim=3 level=3 LinearBoundaryGrid (regular) point examples:
 ...
 161: 1 1 0 | i=0 l=0
 162: 1 1 1 | i=1 l=0
 163: 1 1 2 | i=1 l=0
 164: 1 1 2 | i=1 l=0
 ...
 203: 2 1 2 | i=1 l=1 < coarsable
 204: 2 1 2 | i=1 l=1 < coarsable
 205: 2 2 0 | i=0 l=0
 206: 2 2 1 | i=1 l=1 < coarsable
 207: 2 2 0 | i=0 l=0
 208: 2 2 1 | i=1 l=1 < coarsable
 209: 3 0 0 | i=0 l=0
 210: 3 0 1 | i=0 l=0
 211: 3 1 0 | i=0 l=0
 ...
*/

BOOST_AUTO_TEST_SUITE(TestSurplusVolumeCoarseningFunctor)

/**
   Test basic functionallity
 */
BOOST_AUTO_TEST_CASE(testCoarseningBasic) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  std::vector<size_t>
      toBeRemoved;  // contains grid point hash, not seq numbers!
  DataVector alpha(gridStorage.getSize(), 1.0);

  alpha[161] = 0.1;  // should not be coarsed i=0 l=0
  alpha[162] = 0.1;  // should not be coarsed i=1 l=0
  // should be coarsed i=1 l=1
  alpha[203] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(203).getHash());
  // should be coarsed i=1 l=1
  alpha[204] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(204).getHash());
  // should be coarsed i=1 l=1
  alpha[208] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(208).getHash());

  std::vector<size_t> before;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    before.push_back(it->first->getHash());
  }

  HashCoarsening coarsen;
  SurplusCoarseningFunctor functor(alpha, 3, 0.5);
  coarsen.free_coarsen(gridStorage, functor);

  std::vector<size_t> after;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    after.push_back(it->first->getHash());
  }

  // Check that three points were removed
  BOOST_CHECK_EQUAL(before.size(), after.size() + toBeRemoved.size());

  // Check that non-to-be-removed points are still the same
  for (size_t i = 0; i < before.size(); i++) {
    // check if the point is no longer in the grid storage
    if (std::find(after.begin(), after.end(), before.at(i)) == after.end()) {
      // if so it must have been a removed point
      BOOST_CHECK(std::find(toBeRemoved.begin(), toBeRemoved.end(),
                            before.at(i)) != toBeRemoved.end());
    }
  }

  // Check that all three points really are removed
  for (size_t i = 0; i < toBeRemoved.size(); i++) {
    BOOST_CHECK(std::find(after.begin(), after.end(), toBeRemoved.at(i)) ==
                after.end());
  }
}

/**
 * Test that makes sure the issue with bad_alloc in coarsening does not
 * occur anymore
 */
BOOST_AUTO_TEST_CASE(testBadAllocIssue) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);
  DataVector alpha(gridStorage.getSize(), 0.5);

  HashCoarsening coarsen;
  SurplusCoarseningFunctor functor(alpha, INT_MAX, 1.0);
  coarsen.free_coarsen(gridStorage, functor);
}

/**
 * Test that points with absoulte surplus above the threshold are not coarsed
 */
BOOST_AUTO_TEST_CASE(testCoarseningThreshold) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  std::vector<size_t>
      toBeRemoved;  // contains grid point hash, not seq numbers!
  DataVector alpha(gridStorage.getSize(), 1.0);

  alpha[161] = 0.1;  // should not be coarsed i=0 l=0
  alpha[162] = 0.1;  // should not be coarsed i=1 l=0
  alpha[203] = 0.3;  // should not be coarsed b/c thresh
  // should be coarsed i=1 l=1
  alpha[204] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(204).getHash());
  // should be coarsed i=1 l=1
  alpha[208] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(208).getHash());

  std::vector<size_t> before;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    before.push_back(it->first->getHash());
  }

  HashCoarsening coarsen;
  // threshold of 0.2, max 5 points coarsed
  SurplusCoarseningFunctor functor(alpha, 5, 0.2);
  coarsen.free_coarsen(gridStorage, functor);

  std::vector<size_t> after;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    after.push_back(it->first->getHash());
  }

  // Check that one point was removed
  BOOST_CHECK_EQUAL(before.size(), after.size() + toBeRemoved.size());

  // Check that non-to-be-removed points are still the same
  for (size_t i = 0; i < before.size(); i++) {
    // check if the point is no longer in the grid storage
    if (std::find(after.begin(), after.end(), before.at(i)) == after.end()) {
      // if so it must have been a removed point
      BOOST_CHECK(std::find(toBeRemoved.begin(), toBeRemoved.end(),
                            before.at(i)) != toBeRemoved.end());
    }
  }

  // Check that all three points really are removed
  for (size_t i = 0; i < toBeRemoved.size(); i++) {
    BOOST_CHECK(std::find(after.begin(), after.end(), toBeRemoved.at(i)) ==
                after.end());
  }
}

/**
 * Test that the removedPoints vector filled by the coarsening functionallity
 * really contains
 * exactly the removed points
 * Test further, that removedSeq contains the seq numbers which were removed
 */
BOOST_AUTO_TEST_CASE(testCoarseningRemovedPoints) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  std::vector<size_t>
      toBeRemoved;  // contains grid point hash, not seq numbers!
  std::vector<size_t>
      toBeRemovedSeq;  // contains seq numbers of grid points to be removed
  toBeRemovedSeq.push_back(203);
  toBeRemovedSeq.push_back(204);
  toBeRemovedSeq.push_back(208);
  DataVector alpha(gridStorage.getSize(), 1.0);

  alpha[161] = 0.1;  // should not be coarsed i=0 l=0
  alpha[162] = 0.1;  // should not be coarsed i=1 l=0
  // should be coarsed i=1 l=1
  alpha[203] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(203).getHash());
  // should be coarsed i=1 l=1
  alpha[204] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(204).getHash());
  // should be coarsed i=1 l=1
  alpha[208] = 0.1;
  toBeRemoved.push_back(gridStorage.getPoint(208).getHash());

  std::vector<size_t> before;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    before.push_back(it->first->getHash());
  }

  std::vector<size_t> removedSeq;
  std::vector<HashGridPoint> removedPoints;

  HashCoarsening coarsen;
  SurplusCoarseningFunctor functor(alpha, 3, 0.5);
  coarsen.free_coarsen(gridStorage, functor, &removedPoints, &removedSeq);

  std::vector<size_t> after;
  for (auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    after.push_back(it->first->getHash());
  }
  // Check that toBeRemoved and removedPoints contain the same hash of grid
  // points
  BOOST_CHECK_EQUAL(toBeRemoved.size(), removedPoints.size());
  for (size_t i = 0; i < removedPoints.size(); i++) {
    size_t hash = removedPoints.at(i).getHash();
    bool okay = false;
    for (size_t j = 0; j < toBeRemoved.size(); j++) {
      okay |= (toBeRemoved.at(i) == hash);  // compare the hashs
    }
    BOOST_CHECK(okay);
  }

  // Check that toBeRemovedSeq and removedSeq contain the same elements
  BOOST_CHECK_EQUAL(toBeRemovedSeq.size(), removedSeq.size());
  for (size_t i = 0; i < toBeRemovedSeq.size(); i++) {
    BOOST_CHECK(std::find(removedSeq.begin(), removedSeq.end(),
                          toBeRemovedSeq.at(i)) != removedSeq.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()
