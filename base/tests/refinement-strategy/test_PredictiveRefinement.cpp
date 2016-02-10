// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#include "sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>


using namespace SGPP::base;

BOOST_AUTO_TEST_SUITE(TestPredictiveRefinement)

BOOST_AUTO_TEST_CASE(testFreeRefine2d) {

  size_t dataset_size = 81;
  size_t dim = 2;
  size_t level = 2;

  DataVector error(dataset_size);
  DataMatrix data(dataset_size, dim);

  data.set(0, 0, 0.1);
  data.set(0, 1, 0.1);
  error.set(0, -0.80);
  data.set(1, 0, 0.1);
  data.set(1, 1, 0.2);
  error.set(1, -1.60);
  data.set(2, 0, 0.1);
  data.set(2, 1, 0.3);
  error.set(2, -2.40);
  data.set(3, 0, 0.1);
  data.set(3, 1, 0.4);
  error.set(3, -3.20);
  data.set(4, 0, 0.1);
  data.set(4, 1, 0.5);
  error.set(4, -4.00);
  data.set(5, 0, 0.1);
  data.set(5, 1, 0.6);
  error.set(5, -3.20);
  data.set(6, 0, 0.1);
  data.set(6, 1, 0.7);
  error.set(6, -2.40);
  data.set(7, 0, 0.1);
  data.set(7, 1, 0.8);
  error.set(7, -1.60);
  data.set(8, 0, 0.1);
  data.set(8, 1, 0.9);
  error.set(8, -0.80);
  data.set(9, 0, 0.2);
  data.set(9, 1, 0.1);
  error.set(9, -0.40);
  data.set(10, 0, 0.2);
  data.set(10, 1, 0.2);
  error.set(10, -0.80);
  data.set(11, 0, 0.2);
  data.set(11, 1, 0.3);
  error.set(11, -1.20);
  data.set(12, 0, 0.2);
  data.set(12, 1, 0.4);
  error.set(12, -1.60);
  data.set(13, 0, 0.2);
  data.set(13, 1, 0.5);
  error.set(13, -2.00);
  data.set(14, 0, 0.2);
  data.set(14, 1, 0.6);
  error.set(14, -1.60);
  data.set(15, 0, 0.2);
  data.set(15, 1, 0.7);
  error.set(15, -1.20);
  data.set(16, 0, 0.2);
  data.set(16, 1, 0.8);
  error.set(16, -0.80);
  data.set(17, 0, 0.2);
  data.set(17, 1, 0.9);
  error.set(17, -0.40);
  data.set(18, 0, 0.3);
  data.set(18, 1, 0.1);
  error.set(18, -0.40);
  data.set(19, 0, 0.3);
  data.set(19, 1, 0.2);
  error.set(19, -0.80);
  data.set(20, 0, 0.3);
  data.set(20, 1, 0.3);
  error.set(20, -1.20);
  data.set(21, 0, 0.3);
  data.set(21, 1, 0.4);
  error.set(21, -1.60);
  data.set(22, 0, 0.3);
  data.set(22, 1, 0.5);
  error.set(22, -2.00);
  data.set(23, 0, 0.3);
  data.set(23, 1, 0.6);
  error.set(23, -1.60);
  data.set(24, 0, 0.3);
  data.set(24, 1, 0.7);
  error.set(24, -1.20);
  data.set(25, 0, 0.3);
  data.set(25, 1, 0.8);
  error.set(25, -0.80);
  data.set(26, 0, 0.3);
  data.set(26, 1, 0.9);
  error.set(26, -0.40);
  data.set(27, 0, 0.4);
  data.set(27, 1, 0.1);
  error.set(27, -0.80);
  data.set(28, 0, 0.4);
  data.set(28, 1, 0.2);
  error.set(28, -1.60);
  data.set(29, 0, 0.4);
  data.set(29, 1, 0.3);
  error.set(29, -2.40);
  data.set(30, 0, 0.4);
  data.set(30, 1, 0.4);
  error.set(30, -3.20);
  data.set(31, 0, 0.4);
  data.set(31, 1, 0.5);
  error.set(31, -4.00);
  data.set(32, 0, 0.4);
  data.set(32, 1, 0.6);
  error.set(32, -3.20);
  data.set(33, 0, 0.4);
  data.set(33, 1, 0.7);
  error.set(33, -2.40);
  data.set(34, 0, 0.4);
  data.set(34, 1, 0.8);
  error.set(34, -1.60);
  data.set(35, 0, 0.4);
  data.set(35, 1, 0.9);
  error.set(35, -0.80);
  data.set(36, 0, 0.5);
  data.set(36, 1, 0.1);
  error.set(36, 0.00);
  data.set(37, 0, 0.5);
  data.set(37, 1, 0.2);
  error.set(37, 0.00);
  data.set(38, 0, 0.5);
  data.set(38, 1, 0.3);
  error.set(38, 0.00);
  data.set(39, 0, 0.5);
  data.set(39, 1, 0.4);
  error.set(39, 0.00);
  data.set(40, 0, 0.5);
  data.set(40, 1, 0.5);
  error.set(40, 0.00);
  data.set(41, 0, 0.5);
  data.set(41, 1, 0.6);
  error.set(41, 0.00);
  data.set(42, 0, 0.5);
  data.set(42, 1, 0.7);
  error.set(42, 0.00);
  data.set(43, 0, 0.5);
  data.set(43, 1, 0.8);
  error.set(43, 0.00);
  data.set(44, 0, 0.5);
  data.set(44, 1, 0.9);
  error.set(44, 0.00);
  data.set(45, 0, 0.6);
  data.set(45, 1, 0.1);
  error.set(45, 0.80);
  data.set(46, 0, 0.6);
  data.set(46, 1, 0.2);
  error.set(46, 1.60);
  data.set(47, 0, 0.6);
  data.set(47, 1, 0.3);
  error.set(47, 1.60);
  data.set(48, 0, 0.6);
  data.set(48, 1, 0.4);
  error.set(48, 0.80);
  data.set(49, 0, 0.6);
  data.set(49, 1, 0.5);
  error.set(49, 0.00);
  data.set(50, 0, 0.6);
  data.set(50, 1, 0.6);
  error.set(50, 0.80);
  data.set(51, 0, 0.6);
  data.set(51, 1, 0.7);
  error.set(51, 1.60);
  data.set(52, 0, 0.6);
  data.set(52, 1, 0.8);
  error.set(52, 1.60);
  data.set(53, 0, 0.6);
  data.set(53, 1, 0.9);
  error.set(53, 0.80);
  data.set(54, 0, 0.7);
  data.set(54, 1, 0.1);
  error.set(54, 1.60);
  data.set(55, 0, 0.7);
  data.set(55, 1, 0.2);
  error.set(55, 3.20);
  data.set(56, 0, 0.7);
  data.set(56, 1, 0.3);
  error.set(56, 3.20);
  data.set(57, 0, 0.7);
  data.set(57, 1, 0.4);
  error.set(57, 1.60);
  data.set(58, 0, 0.7);
  data.set(58, 1, 0.5);
  error.set(58, 0.00);
  data.set(59, 0, 0.7);
  data.set(59, 1, 0.6);
  error.set(59, 1.60);
  data.set(60, 0, 0.7);
  data.set(60, 1, 0.7);
  error.set(60, 3.20);
  data.set(61, 0, 0.7);
  data.set(61, 1, 0.8);
  error.set(61, 3.20);
  data.set(62, 0, 0.7);
  data.set(62, 1, 0.9);
  error.set(62, 1.60);
  data.set(63, 0, 0.8);
  data.set(63, 1, 0.1);
  error.set(63, 1.60);
  data.set(64, 0, 0.8);
  data.set(64, 1, 0.2);
  error.set(64, 3.20);
  data.set(65, 0, 0.8);
  data.set(65, 1, 0.3);
  error.set(65, 3.20);
  data.set(66, 0, 0.8);
  data.set(66, 1, 0.4);
  error.set(66, 1.60);
  data.set(67, 0, 0.8);
  data.set(67, 1, 0.5);
  error.set(67, 0.00);
  data.set(68, 0, 0.8);
  data.set(68, 1, 0.6);
  error.set(68, 1.60);
  data.set(69, 0, 0.8);
  data.set(69, 1, 0.7);
  error.set(69, 3.20);
  data.set(70, 0, 0.8);
  data.set(70, 1, 0.8);
  error.set(70, 3.20);
  data.set(71, 0, 0.8);
  data.set(71, 1, 0.9);
  error.set(71, 1.60);
  data.set(72, 0, 0.9);
  data.set(72, 1, 0.1);
  error.set(72, 0.80);
  data.set(73, 0, 0.9);
  data.set(73, 1, 0.2);
  error.set(73, 1.60);
  data.set(74, 0, 0.9);
  data.set(74, 1, 0.3);
  error.set(74, 1.60);
  data.set(75, 0, 0.9);
  data.set(75, 1, 0.4);
  error.set(75, 0.80);
  data.set(76, 0, 0.9);
  data.set(76, 1, 0.5);
  error.set(76, 0.00);
  data.set(77, 0, 0.9);
  data.set(77, 1, 0.6);
  error.set(77, 0.80);
  data.set(78, 0, 0.9);
  data.set(78, 1, 0.7);
  error.set(78, 1.60);
  data.set(79, 0, 0.9);
  data.set(79, 1, 0.8);
  error.set(79, 1.60);
  data.set(80, 0, 0.9);
  data.set(80, 1, 0.9);
  error.set(80, 0.80);


  Grid* grid = Grid::createLinearGrid(dim);

  GridGenerator* gen = grid->createGridGenerator();
  gen->regular(level);
  //  GridStorage* storage = grid->getStorage();

  DataVector surplusses(1);
  surplusses[0] = 1.0;
  surplusses[1] = .5;
  surplusses[2] = .5;
  surplusses[3] = .5;
  surplusses[4] = .5;

  //  PredictiveRefinementIndicator functor(grid, &data, &error, 2, 0., 0. );
  //  HashRefinement* hash_refinement = new HashRefinement();
  //  PredictiveRefinement predictive_refinement(hash_refinement);
  //
  //
  //  predictive_refinement.free_refine(storage, &functor);
  //
  //  BOOST_CHECK_EQUAL(storage->size(), 9);
  //
  //  for (size_t i = 0; i < storage->size(); i++) {
  //    HashGridIndex* index = storage->get(i);
  //    BOOST_CHECK((index->getIndex(0) == 4) == false);
  //  }
  //
  //
  //  delete hash_refinement;
  //  delete gen;
  //  delete grid;
}


BOOST_AUTO_TEST_SUITE_END()
