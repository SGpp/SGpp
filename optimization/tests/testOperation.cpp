
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/optimization/function/test/Sphere.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include "GridCreator.hpp"

const bool use_double_precision =
#if USE_DOUBLE_PRECISION
  true;
#else
  false;
#endif /* USE_DOUBLE_PRECISION */

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestOperationMultipleHierarchisation) {
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 3;
  const size_t p = 5;
  const size_t l = 4;
  const size_t m = 4;
  const SGPP::float_t tol = (use_double_precision ? 1e-4 : 5e-3);

  test_functions::Sphere f(d);

  // Test All The Grids!
  std::vector<std::unique_ptr<base::Grid>> grids;
  createSupportedGrids(d, p, grids);

  for (auto& grid : grids) {
    base::DataVector functionValues(0);
    f.generateDisplacement();
    createSampleGrid(*grid, l, f, functionValues);

    base::DataMatrix functionValuesMatrix(grid->getSize(), m);

    for (size_t j = 0; j < m; j++) {
      base::DataVector column(0);
      f.generateDisplacement();
      createSampleGrid(*grid, l, f, column);
      functionValuesMatrix.setColumn(j, column);
    }

    std::unique_ptr<OperationMultipleHierarchisation> op(
      op_factory::createOperationMultipleHierarchisation(*grid));

    base::DataVector functionValues2(functionValues);
    op->doHierarchisation(functionValues2);
    op->doDehierarchisation(functionValues2);

    for (size_t i = 0; i < grid->getSize(); i++) {
      BOOST_CHECK_CLOSE(functionValues[i], functionValues2[i], tol);
    }

    base::DataMatrix functionValuesMatrix2(functionValuesMatrix);
    op->doHierarchisation(functionValuesMatrix2);
    op->doDehierarchisation(functionValuesMatrix2);

    for (size_t i = 0; i < grid->getSize(); i++) {
      for (size_t j = 0; j < m; j++) {
        BOOST_CHECK_CLOSE(functionValuesMatrix.get(i, j),
                          functionValuesMatrix2.get(i, j), tol);
      }
    }
  }
}
