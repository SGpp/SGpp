// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp>

#include <iostream>
#include <algorithm>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::HashRefinementBoundaries;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::SubspaceRefinement;
using sgpp::base::AbstractRefinement;
using sgpp::base::RefinementFunctor;


BOOST_AUTO_TEST_SUITE(TestRefinementAddedPoints)

/*
  Utility for tests
 */
void gridTest(GridStorage& gridStorage, DataVector& alpha, AbstractRefinement& ref) {
  std::vector<size_t> oldPoints; // points before refinement
  std::vector<size_t> newPoints; // points that are new after refinement
  std::vector<size_t> addedPoints; // points that are added by refinement
  for(auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    oldPoints.push_back(gridStorage.getSequenceNumber(*it->first));
  }

  SurplusRefinementFunctor fun(alpha, 2);

  ref.free_refine(gridStorage, fun, &addedPoints);

  for(auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    size_t seq = gridStorage.getSequenceNumber(*it->first);
    // only added points not contained in oldPoints to newPoints
    if(std::find(oldPoints.begin(), oldPoints.end(),seq) == oldPoints.end()) {
      newPoints.push_back(seq);
    }
  }

  // check if all newPoints == addedPoints (up to order)
  BOOST_CHECK_EQUAL(addedPoints.size(), newPoints.size());
  for(size_t i = 0; i < addedPoints.size(); i++) {
    bool b = std::find(newPoints.begin(), newPoints.end(), addedPoints[i])
      != newPoints.end();
    BOOST_CHECK(b);
  }
}

BOOST_AUTO_TEST_CASE(TestBasicHash) {
  size_t dim = 4;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  // arbitrary choice
  alpha[3] = 1.0; // level 3 1 1 1
  alpha[37] = 1.0; // level 2 1 1 2

  HashRefinement refHash;
  gridTest(gridStorage, alpha, refHash);
}

BOOST_AUTO_TEST_CASE(TestBasicSubspace) {
  size_t dim = 4;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  // arbitrary choice
  alpha[3] = 1.0; // level 3 1 1 1
  alpha[37] = 1.0; // level 2 1 1 2

  HashRefinement refHash;
  SubspaceRefinement refSub(&refHash);
  gridTest(gridStorage, alpha, refSub);
}

BOOST_AUTO_TEST_CASE(TestLinearPredictive) {
  // TODO
}

BOOST_AUTO_TEST_CASE(TestBoundaryHash) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  // arbitrary choice
  alpha[30] = 1.0; // level 1 3 0
  alpha[180] = 1.0; // level 1 3 1

  HashRefinementBoundaries refHash;
  gridTest(gridStorage, alpha, refHash);
}

BOOST_AUTO_TEST_CASE(TestBoundarySubspace) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  // arbitrary choice
  alpha[30] = 1.0; // level 1 3 0
  alpha[180] = 1.0; // level 1 3 1

  HashRefinementBoundaries refHash;
  SubspaceRefinement refSub(&refHash);
  gridTest(gridStorage, alpha, refSub);
}

BOOST_AUTO_TEST_CASE(TestBoundaryPredictive) {
  // TODO
}

BOOST_AUTO_TEST_SUITE_END()
