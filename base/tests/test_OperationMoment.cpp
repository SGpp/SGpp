// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

double epsilon = 1e-6;        // maximum difference between the approximation and result

double firstMomentApproximation(base::Grid* grid, base::DataVector& alpha) {
  const size_t resolution = 10000;
  const double h = 1.0 / static_cast<double>(resolution);
  base::GridStorage& storage = grid->getStorage();
  if (storage.getSize() != alpha.getSize()) {
    throw base::application_exception("FirstMomentApproximation: Grid and alpha don't match");
  }
  base::SBasis& basis = const_cast<sgpp::base::SBasis&>(grid->getBasis());
  // trapezoidal rule

  double res = 0.0;
  for (size_t i = 0; i < grid->getSize(); i++) {
    const sgpp::base::level_t level = storage.getPoint(i).getLevel(0);
    const sgpp::base::index_t index = storage.getPoint(i).getIndex(0);
    double temp_res = 0.0;
    // --------------------------------------------------------------------------
    // apply trapezoidal rule
    // x = 0 term can be neglected
    for (size_t c = 1; c < resolution; c++) {
      double x = static_cast<double>(c) * h;
      temp_res += x * basis.eval(level, index, x);
    }
    // x = 1. term
    temp_res += basis.eval(level, index, 0.0) / 2.0;
    temp_res *= h;
    // --------------------------------------------------------------------------
    res += alpha.get(i) * temp_res;
  }
  return res;
}

BOOST_AUTO_TEST_SUITE(testOperationFirstMoment)

BOOST_AUTO_TEST_CASE(testOperationFirstMomentLinear) {
  const size_t d = 1;
  const size_t l = 5;
  base::Grid* grid(sgpp::base::Grid::createLinearGrid(d));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentModLinear) {
  const size_t d = 1;
  const size_t l = 5;
  base::Grid* grid(sgpp::base::Grid::createModLinearGrid(d));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentPoly) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createPolyGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentModPoly) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createModPolyGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentPolyBoundary) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createPolyBoundaryGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentPolyClenshwCurtis) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentModPolyClenshawCurtis) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createModPolyClenshawCurtisGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentPolyClenshawCurtisBoundary) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisBoundaryGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}


BOOST_AUTO_TEST_CASE(testOperationFirstMomentBspline) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createBsplineGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentModBspline) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createModBsplineGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentBsplineBoundary) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createBsplineBoundaryGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentBsplineClenshwCurtis) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_CASE(testOperationFirstMomentModBsplineClenshawCurtis) {
  const size_t d = 1;
  const size_t l = 4;
  base::Grid* grid(sgpp::base::Grid::createModBsplineClenshawCurtisGrid(d, 3));
  grid->getGenerator().regular(l);
  base::OperationFirstMoment* op = sgpp::op_factory::createOperationFirstMoment(*grid);
  base::DataVector alpha(grid->getSize(), 1.0);
  double approx = firstMomentApproximation(grid, alpha);
  double res = op->doQuadrature(alpha, nullptr);
  BOOST_CHECK_SMALL(approx - res, epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace pde
}  // namespace sgpp
