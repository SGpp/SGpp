// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/MultipleClassRefinement.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>

#include <vector>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;

BOOST_AUTO_TEST_SUITE(MultipleClassRefinementFunctor)

BOOST_AUTO_TEST_CASE(testRefine2class) {
    DataVector alphas1(50);
    alphas1.set(0, 0.05);
    alphas1.set(1, -0.05);
    alphas1.set(2, 0.05);
    alphas1.set(3, -0.05);
    alphas1.set(4, 0.05);
    alphas1.set(5, -0.05);
    alphas1.set(6, 0.05);
    alphas1.set(7, -0.05);
    alphas1.set(8, 0.05);
    alphas1.set(9, -0.05);
    alphas1.set(10, -0.05);
    alphas1.set(11, -0.05);
    alphas1.set(12, 0.05);
    alphas1.set(13, 0.05);
    alphas1.set(14, -0.05);
    alphas1.set(15, 0.05);
    alphas1.set(16, 0.05);
    alphas1.set(17, 0.05);
    alphas1.set(18, -0.05);
    alphas1.set(19, -0.05);
    alphas1.set(20, 0.05);
    alphas1.set(21, -0.05);
    alphas1.set(22, 0.05);
    alphas1.set(23, 0.05);
    alphas1.set(24, -0.05);
    alphas1.set(25, -0.05);
    alphas1.set(26, -0.05);
    alphas1.set(27, 0.05);
    alphas1.set(28, -0.05);
    alphas1.set(29, -0.05);
    alphas1.set(30, 0.05);
    alphas1.set(31, -0.05);
    alphas1.set(32, 0.05);
    alphas1.set(33, 0.05);
    alphas1.set(34, 0.05);
    alphas1.set(35, 0.05);
    alphas1.set(36, -0.05);
    alphas1.set(37, -0.05);
    alphas1.set(38, -0.05);
    alphas1.set(39, -0.05);
    alphas1.set(40, 0.05);
    alphas1.set(41, 0.05);
    alphas1.set(42, -0.05);
    alphas1.set(43, -0.05);
    alphas1.set(44, -0.05);
    alphas1.set(45, 0.05);
    alphas1.set(46, 0.05);
    alphas1.set(47, -0.05);
    alphas1.set(48, -0.05);
    alphas1.set(49, 0.05);

    DataVector alphas2(50);
    alphas2.set(0, -0.05);
    alphas2.set(1, 0.05);
    alphas2.set(2, -0.05);
    alphas2.set(3, 0.05);
    alphas2.set(4, -0.05);
    alphas2.set(5, 0.05);
    alphas2.set(6, -0.05);
    alphas2.set(7, 0.05);
    alphas2.set(8, -0.05);
    alphas2.set(9, 0.05);
    alphas2.set(10, 0.05);
    alphas2.set(11, 0.05);
    alphas2.set(12, -0.05);
    alphas2.set(13, -0.05);
    alphas2.set(14, 0.05);
    alphas2.set(15, -0.05);
    alphas2.set(16, -0.05);
    alphas2.set(17, -0.05);
    alphas2.set(18, 0.05);
    alphas2.set(19, 0.05);
    alphas2.set(20, -0.05);
    alphas2.set(21, 0.05);
    alphas2.set(22, -0.05);
    alphas2.set(23, -0.05);
    alphas2.set(24, 0.05);
    alphas2.set(25, 0.05);
    alphas2.set(26, 0.05);
    alphas2.set(27, -0.05);
    alphas2.set(28, 0.05);
    alphas2.set(29, 0.05);
    alphas2.set(30, -0.05);
    alphas2.set(31, 0.05);
    alphas2.set(32, -0.05);
    alphas2.set(33, -0.05);
    alphas2.set(34, -0.05);
    alphas2.set(35, -0.05);
    alphas2.set(36, 0.05);
    alphas2.set(37, 0.05);
    alphas2.set(38, 0.05);
    alphas2.set(39, 0.05);
    alphas2.set(40, -0.05);
    alphas2.set(41, -0.05);
    alphas2.set(42, 0.05);
    alphas2.set(43, 0.05);
    alphas2.set(44, 0.05);
    alphas2.set(45, -0.05);
    alphas2.set(46, -0.05);
    alphas2.set(47, 0.05);
    alphas2.set(48, 0.05);
    alphas2.set(49, -0.05);

    size_t dim = 2;
    size_t level = 3;

    std::vector<sgpp::base::Grid*> grids;
    std::vector<sgpp::base::DataVector*> alphas;
    std::vector<double> priors;

    std::shared_ptr<sgpp::base::Grid> grid1(sgpp::base::Grid::createLinearGrid(dim));
    grid1->getGenerator().regular(level);
    grids.push_back(grid1.get());
    alphas.push_back(&alphas1);
    priors.push_back(1.0);

    std::shared_ptr<sgpp::base::Grid> grid2(sgpp::base::Grid::createLinearGrid(dim));
    grid2->getGenerator().regular(level);
    grids.push_back(grid2.get());
    alphas.push_back(&alphas2);
    priors.push_back(1.0);

    size_t sizeGrids = 17;

    BOOST_CHECK_EQUAL(grid1->getStorage().getSize(), sizeGrids);
    BOOST_CHECK_EQUAL(grid2->getStorage().getSize(), sizeGrids);

    sgpp::datadriven::MultipleClassRefinementFunctor mcrf(grids, alphas, priors, 2, 0, 0);
    mcrf.refine();

    BOOST_CHECK_GT(grid1->getStorage().getSize(), sizeGrids);
    BOOST_CHECK_GT(grid2->getStorage().getSize(), sizeGrids);
}

BOOST_AUTO_TEST_SUITE_END()
