// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::HashGridPoint;
using sgpp::base::HashRefinement;
using sgpp::base::ImpurityRefinement;
using sgpp::base::ImpurityRefinementIndicator;

BOOST_AUTO_TEST_SUITE(TestImpurityRefinement)

BOOST_AUTO_TEST_CASE(testFreeRefine2d) {
  DataMatrix data(50, 2);
  data.set(0, 0, 0.561286);
  data.set(0, 1, 0.37669);
  data.set(1, 0, 0.831794);
  data.set(1, 1, 0.347549);
  data.set(2, 0, 0.785171);
  data.set(2, 1, 0.481398);
  data.set(3, 0, 0.374007);
  data.set(3, 1, 0.681501);
  data.set(4, 0, 0.00495);
  data.set(4, 1, 0.476371);
  data.set(5, 0, 0.458995);
  data.set(5, 1, 0.785822);
  data.set(6, 0, 0.581328);
  data.set(6, 1, 0.379036);
  data.set(7, 0, 0.442863);
  data.set(7, 1, 0.654657);
  data.set(8, 0, 0.633425);
  data.set(8, 1, 0.470952);
  data.set(9, 0, 0.560704);
  data.set(9, 1, 0.434832);
  data.set(10, 0, 0.663335);
  data.set(10, 1, 0.426796);
  data.set(11, 0, 0.835763);
  data.set(11, 1, 0.448754);
  data.set(12, 0, 0.690902);
  data.set(12, 1, 0.394863);
  data.set(13, 0, 0.717678);
  data.set(13, 1, 0.334948);
  data.set(14, 0, 0.707282);
  data.set(14, 1, 0.656742);
  data.set(15, 0, 0.660405);
  data.set(15, 1, 0.664967);
  data.set(16, 0, 0.293661);
  data.set(16, 1, 0.362652);
  data.set(17, 0, 0.826191);
  data.set(17, 1, 0.217398);
  data.set(18, 0, 0.493762);
  data.set(18, 1, 0.659877);
  data.set(19, 0, 0.466227);
  data.set(19, 1, 0.53704);
  data.set(20, 0, 0.688684);
  data.set(20, 1, 0.487203);
  data.set(21, 0, 0.725357);
  data.set(21, 1, 0.469516);
  data.set(22, 0, 0.36625);
  data.set(22, 1, 0.612607);
  data.set(23, 0, 0.583864);
  data.set(23, 1, 0.582546);
  data.set(24, 0, 0.352386);
  data.set(24, 1, 0.70859);
  data.set(25, 0, 0.723137);
  data.set(25, 1, 0.691439);
  data.set(26, 0, 0.672744);
  data.set(26, 1, 0.466487);
  data.set(27, 0, 0.696203);
  data.set(27, 1, 0.64423);
  data.set(28, 0, 0.811522);
  data.set(28, 1, 0.696915);
  data.set(29, 0, 0.521949);
  data.set(29, 1, 0.548102);
  data.set(30, 0, 0.358435);
  data.set(30, 1, 0.310331);
  data.set(31, 0, 0.724029);
  data.set(31, 1, 0.61844);
  data.set(32, 0, 0.747972);
  data.set(32, 1, 0.443816);
  data.set(33, 0, 0.616983);
  data.set(33, 1, 0.27329);
  data.set(34, 0, 0.658415);
  data.set(34, 1, 0.517102);
  data.set(35, 0, 0.601694);
  data.set(35, 1, 0.563502);
  data.set(36, 0, 0.727661);
  data.set(36, 1, 0.708638);
  data.set(37, 0, 0.788582);
  data.set(37, 1, 0.537601);
  data.set(38, 0, 0.749369);
  data.set(38, 1, 0.564056);
  data.set(39, 0, 0.41065);
  data.set(39, 1, 0.432474);
  data.set(40, 0, 0.746131);
  data.set(40, 1, 0.23194);
  data.set(41, 0, 0.811588);
  data.set(41, 1, 0.485242);
  data.set(42, 0, 0.454642);
  data.set(42, 1, 0.721785);
  data.set(43, 0, 0.781686);
  data.set(43, 1, 0.593543);
  data.set(44, 0, 0.718447);
  data.set(44, 1, 0.504362);
  data.set(45, 0, 0.17071);
  data.set(45, 1, 0.338819);
  data.set(46, 0, 0.576735);
  data.set(46, 1, 0.260155);
  data.set(47, 0, 0.758747);
  data.set(47, 1, 0.513633);
  data.set(48, 0, 0.695441);
  data.set(48, 1, 0.531488);
  data.set(49, 0, 0.335224);
  data.set(49, 1, 0.578003);

  DataVector classes(50);
  classes.set(0, -1);
  classes.set(1, -1);
  classes.set(2, 1);
  classes.set(3, 1);
  classes.set(4, -1);
  classes.set(5, 1);
  classes.set(6, -1);
  classes.set(7, 1);
  classes.set(8, -1);
  classes.set(9, -1);
  classes.set(10, -1);
  classes.set(11, 1);
  classes.set(12, -1);
  classes.set(13, -1);
  classes.set(14, 1);
  classes.set(15, 1);
  classes.set(16, -1);
  classes.set(17, -1);
  classes.set(18, 1);
  classes.set(19, 1);
  classes.set(20, 1);
  classes.set(21, 1);
  classes.set(22, 1);
  classes.set(23, 1);
  classes.set(24, 1);
  classes.set(25, 1);
  classes.set(26, -1);
  classes.set(27, 1);
  classes.set(28, 1);
  classes.set(29, 1);
  classes.set(30, -1);
  classes.set(31, 1);
  classes.set(32, 1);
  classes.set(33, -1);
  classes.set(34, 1);
  classes.set(35, 1);
  classes.set(36, 1);
  classes.set(37, 1);
  classes.set(38, 1);
  classes.set(39, -1);
  classes.set(40, -1);
  classes.set(41, 1);
  classes.set(42, 1);
  classes.set(43, 1);
  classes.set(44, 1);
  classes.set(45, -1);
  classes.set(46, -1);
  classes.set(47, 1);
  classes.set(48, 1);
  classes.set(49, 1);

  DataMatrix gridPoints(27, 4);
  gridPoints.set(0, 0, 3);
  gridPoints.set(0, 1, 2);
  gridPoints.set(0, 2, 5);
  gridPoints.set(0, 3, 3);
  gridPoints.set(1, 0, 4);
  gridPoints.set(1, 1, 1);
  gridPoints.set(1, 2, 11);
  gridPoints.set(1, 3, 1);
  gridPoints.set(2, 0, 4);
  gridPoints.set(2, 1, 1);
  gridPoints.set(2, 2, 9);
  gridPoints.set(2, 3, 1);
  gridPoints.set(3, 0, 2);
  gridPoints.set(3, 1, 3);
  gridPoints.set(3, 2, 3);
  gridPoints.set(3, 3, 1);
  gridPoints.set(4, 0, 3);
  gridPoints.set(4, 1, 2);
  gridPoints.set(4, 2, 7);
  gridPoints.set(4, 3, 1);
  gridPoints.set(5, 0, 1);
  gridPoints.set(5, 1, 1);
  gridPoints.set(5, 2, 1);
  gridPoints.set(5, 3, 1);
  gridPoints.set(6, 0, 1);
  gridPoints.set(6, 1, 3);
  gridPoints.set(6, 2, 1);
  gridPoints.set(6, 3, 1);
  gridPoints.set(7, 0, 2);
  gridPoints.set(7, 1, 1);
  gridPoints.set(7, 2, 1);
  gridPoints.set(7, 3, 1);
  gridPoints.set(8, 0, 1);
  gridPoints.set(8, 1, 2);
  gridPoints.set(8, 2, 1);
  gridPoints.set(8, 3, 3);
  gridPoints.set(9, 0, 1);
  gridPoints.set(9, 1, 2);
  gridPoints.set(9, 2, 1);
  gridPoints.set(9, 3, 1);
  gridPoints.set(10, 0, 1);
  gridPoints.set(10, 1, 3);
  gridPoints.set(10, 2, 1);
  gridPoints.set(10, 3, 7);
  gridPoints.set(11, 0, 2);
  gridPoints.set(11, 1, 1);
  gridPoints.set(11, 2, 3);
  gridPoints.set(11, 3, 1);
  gridPoints.set(12, 0, 1);
  gridPoints.set(12, 1, 3);
  gridPoints.set(12, 2, 1);
  gridPoints.set(12, 3, 3);
  gridPoints.set(13, 0, 2);
  gridPoints.set(13, 1, 2);
  gridPoints.set(13, 2, 1);
  gridPoints.set(13, 3, 1);
  gridPoints.set(14, 0, 1);
  gridPoints.set(14, 1, 3);
  gridPoints.set(14, 2, 1);
  gridPoints.set(14, 3, 5);
  gridPoints.set(15, 0, 2);
  gridPoints.set(15, 2, 1);
  gridPoints.set(15, 1, 2);
  gridPoints.set(15, 3, 3);
  gridPoints.set(16, 0, 2);
  gridPoints.set(16, 1, 2);
  gridPoints.set(16, 2, 3);
  gridPoints.set(16, 3, 1);
  gridPoints.set(17, 0, 2);
  gridPoints.set(17, 1, 3);
  gridPoints.set(17, 2, 1);
  gridPoints.set(17, 3, 3);
  gridPoints.set(18, 0, 2);
  gridPoints.set(18, 1, 2);
  gridPoints.set(18, 2, 3);
  gridPoints.set(18, 3, 3);
  gridPoints.set(19, 0, 3);
  gridPoints.set(19, 1, 1);
  gridPoints.set(19, 2, 1);
  gridPoints.set(19, 3, 1);
  gridPoints.set(20, 0, 1);
  gridPoints.set(20, 1, 4);
  gridPoints.set(20, 2, 1);
  gridPoints.set(20, 3, 5);
  gridPoints.set(21, 0, 3);
  gridPoints.set(21, 1, 1);
  gridPoints.set(21, 2, 3);
  gridPoints.set(21, 3, 1);
  gridPoints.set(22, 0, 3);
  gridPoints.set(22, 1, 1);
  gridPoints.set(22, 2, 5);
  gridPoints.set(22, 3, 1);
  gridPoints.set(23, 0, 3);
  gridPoints.set(23, 1, 1);
  gridPoints.set(23, 2, 7);
  gridPoints.set(23, 3, 1);
  gridPoints.set(24, 0, 2);
  gridPoints.set(24, 1, 3);
  gridPoints.set(24, 2, 3);
  gridPoints.set(24, 3, 3);
  gridPoints.set(25, 0, 1);
  gridPoints.set(25, 2, 1);
  gridPoints.set(25, 1, 4);
  gridPoints.set(25, 3, 7);
  gridPoints.set(26, 0, 3);
  gridPoints.set(26, 1, 2);
  gridPoints.set(26, 2, 5);
  gridPoints.set(26, 3, 1);

  std::unique_ptr<Grid> grid;
  grid.reset(Grid::createModLinearGrid(2));
  grid->getGenerator().regular(3);
  GridStorage& gridStorage = grid->getStorage();

  HashRefinement refinement;
  ImpurityRefinement decorator(&refinement);
  ImpurityRefinementIndicator indicator(*grid, data, nullptr, nullptr, nullptr,
                                        classes, 0.0, 3);

  decorator.free_refine(gridStorage, indicator);

  BOOST_CHECK_EQUAL(gridStorage.getSize(), 27);

  HashGridPoint *point = new HashGridPoint(2);
  for (size_t i = 0; i < 27; i++) {
      point->set(0, static_cast<unsigned int>(gridPoints.get(i, 0)),
          static_cast<unsigned int>(gridPoints.get(i, 2)));
      point->set(1, static_cast<unsigned int>(gridPoints.get(i, 1)),
          static_cast<unsigned int>(gridPoints.get(i, 3)));
      BOOST_CHECK(gridStorage.isContaining(*point));
  }
  delete point;
}

BOOST_AUTO_TEST_SUITE_END()
