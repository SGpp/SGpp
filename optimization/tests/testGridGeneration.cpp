// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rosenbrock.hpp>

#include <vector>

#include "GridCreator.hpp"

using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;
using sgpp::base::ScalarFunction;
using sgpp::optimization::IterativeGridGenerator;
using sgpp::optimization::IterativeGridGeneratorLinearSurplus;
using sgpp::optimization::IterativeGridGeneratorRitterNovak;
using sgpp::optimization::IterativeGridGeneratorSOO;
using sgpp::optimization::test_problems::Rosenbrock;

BOOST_AUTO_TEST_CASE(TestIterativeGridGenerators) {
  // Test sgpp::optimization iterative grid generators.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t N = 200;

  Rosenbrock testProblem(d);
  testProblem.generateDisplacement();
  ScalarFunction& f = testProblem.getObjectiveFunction();

  // Test All The Grids!
  std::vector<std::unique_ptr<sgpp::base::Grid>> grids;
  createSupportedGrids(d, p, grids);

  // test getters/setters
  {
    IterativeGridGeneratorRitterNovak gridGen(f, *grids[0], N);

    BOOST_CHECK_EQUAL(&gridGen.getGrid(), grids[0].get());

    const double adaptivity = 0.42;
    gridGen.setAdaptivity(adaptivity);
    BOOST_CHECK_EQUAL(gridGen.getAdaptivity(), adaptivity);

    const sgpp::base::level_t maxLevel = 13;
    gridGen.setMaxLevel(maxLevel);
    BOOST_CHECK_EQUAL(gridGen.getMaxLevel(), maxLevel);

    const IterativeGridGeneratorRitterNovak::PowMethod powMethod =
        IterativeGridGeneratorRitterNovak::PowMethod::FAST_POW;
    gridGen.setPowMethod(powMethod);
    BOOST_CHECK(gridGen.getPowMethod() == powMethod);
  }

  {
    IterativeGridGeneratorLinearSurplus gridGen(f, *grids[0], N);

    BOOST_CHECK_EQUAL(&gridGen.getGrid(), grids[0].get());

    const double adaptivity = 0.42;
    gridGen.setAdaptivity(adaptivity);
    BOOST_CHECK_EQUAL(gridGen.getAdaptivity(), adaptivity);
  }

  {
    IterativeGridGeneratorSOO gridGen(f, *grids[0], N);

    BOOST_CHECK_EQUAL(&gridGen.getGrid(), grids[0].get());

    const double adaptivity = 0.42;
    gridGen.setAdaptivity(adaptivity);

    const IterativeGridGeneratorSOO::AdaptivityFunction adaptivityFunction = [](size_t n) {
      return n * n;
    };
    gridGen.setAdaptivity(adaptivityFunction);
    BOOST_CHECK_EQUAL(gridGen.getAdaptivity()(42), static_cast<size_t>(42 * 42));
  }

  for (auto& grid : grids) {
    // repeat for grid generators
    IterativeGridGeneratorRitterNovak gridGenRN(f, *grid, N, 0.85);
    IterativeGridGeneratorRitterNovak gridGenRNFastPow(f, *grid, N, 0.85);
    IterativeGridGeneratorLinearSurplus gridGenLS(f, *grid, N, 0.85);
    IterativeGridGeneratorSOO gridGenSOO(f, *grid, N, 0.85);

    gridGenRNFastPow.setPowMethod(IterativeGridGeneratorRitterNovak::PowMethod::FAST_POW);

    std::vector<IterativeGridGenerator*> gridGens = {&gridGenRN, &gridGenRNFastPow, &gridGenLS,
                                                     &gridGenSOO};

    for (auto& gridGen : gridGens) {
      // empty grid
      grid->getStorage().clear();

      // generate grid
      BOOST_CHECK(gridGen->generate());

      // test grid size
      const size_t n = grid->getSize();
      BOOST_CHECK_LE(n, N);

      // test size of function value vector
      const sgpp::base::DataVector& functionValues = gridGen->getFunctionValues();
      BOOST_CHECK_EQUAL(n, functionValues.getSize());

      for (size_t i = 0; i < n; i++) {
        sgpp::base::DataVector x(grid->getStorage().getCoordinates((grid->getStorage())[i]));
        // test function value
        BOOST_CHECK_CLOSE(functionValues[i], f.eval(x), 1e-10);
      }

      // test printing of grid points and function values
      gridGen->printIterativeGridGenerator();
    }
  }
}
