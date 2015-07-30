#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
//#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

using namespace SGPP::base;
BOOST_AUTO_TEST_SUITE(TestOperationMultipleEval)

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  size_t dim = 2;
  Grid* grid = Grid::createLinearGrid(dim);
  grid->createGridGenerator()->regular(2);

  GridStorage* gS = grid->getStorage();

  size_t N = gS->size();

  DataVector alpha(N);

  for (int i = 0; i < static_cast<int>(N); ++i) {
    alpha[i] = i + 1;
  }

  SGPP::float_t points[3][2] = { { 0.5, 1.0 }, { 0.3, 0.4 }, { 0.9, 0.7 } };
  size_t numberDataPoints = 3;

  DataVector result(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result[i] = i;
  }

  DataMatrix dataset(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      temp[j] = points[i][j];
    }

    dataset.setRow(i, temp);
  }

  OperationMultipleEval* multiEvalOp =
    SGPP::op_factory::createOperationMultipleEval(*grid, dataset);

  multiEvalOp->mult(alpha, result);

  BOOST_MESSAGE(alpha.toString() + "\n");
  BOOST_MESSAGE(result.toString() + "\n");

  DataVector result_ref(numberDataPoints);
  result_ref[0] = 0.0;
  result_ref[1] = 2.72;
  result_ref[2] = 1.64;

  BOOST_CHECK_CLOSE(result[0], result_ref[0], 1e-7);
  BOOST_CHECK_CLOSE(result[1], result_ref[1], 1e-7);
  BOOST_CHECK_CLOSE(result[2], result_ref[2], 1e-7);

  delete grid;
  delete multiEvalOp;
}

BOOST_AUTO_TEST_SUITE_END()
