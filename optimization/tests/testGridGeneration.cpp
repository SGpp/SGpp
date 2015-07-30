#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/test/Rosenbrock.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestIterativeGridGenerators) {
  // Test SGPP::optimization iterative grid generators.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t N = 200;

  test_functions::Rosenbrock f(d);
  f.generateDisplacement();

  // Test All The Grids!
  std::vector<std::unique_ptr<base::Grid>> grids;
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineTruncatedBoundaryGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createBsplineClenshawCurtisGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModBsplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createFundamentalSplineGrid(d, p))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearTruncatedBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createLinearClenshawCurtisGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModLinearGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createWaveletGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createWaveletTruncatedBoundaryGrid(d))));
  grids.push_back(std::move(std::unique_ptr<base::Grid>(
                              base::Grid::createModWaveletGrid(d))));

  for (auto& grid : grids) {
    // repeat for Ritter-Novak and linear surplus grid generation
    std::vector<std::unique_ptr<IterativeGridGenerator>> gridGens;
    gridGens.push_back(
      std::move(std::unique_ptr<IterativeGridGenerator>(
                  new IterativeGridGeneratorRitterNovak(f, *grid, N, 0.85))));
    gridGens.push_back(
      std::move(std::unique_ptr<IterativeGridGenerator>(
                  new IterativeGridGeneratorLinearSurplus(f, *grid, N, 0.2))));

    for (auto& gridGen : gridGens) {
      // empty grid
      grid->getStorage()->emptyStorage();

      // generate grid
      BOOST_CHECK(gridGen->generate());

      // test grid size
      const size_t n = grid->getSize();
      BOOST_CHECK_LE(n, N);

      // test size of function value vector
      const base::DataVector& functionValues = gridGen->getFunctionValues();
      BOOST_CHECK_EQUAL(n, functionValues.getSize());

      for (size_t i = 0; i < n; i++) {
        base::DataVector x(d);

        for (size_t t = 0; t < d; t++) {
          x[t] = grid->getStorage()->get(i)->getCoord(t);
        }

        // test function value
        BOOST_CHECK_CLOSE(functionValues.get(i), f.eval(x), 1e-10);
      }
    }
  }
}
