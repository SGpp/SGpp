// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>

#include <vector>

using sgpp::base::BoundingBox1D;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::HashGridPoint;
using sgpp::base::OperationQuadrature;
using sgpp::base::OperationQuadratureMC;
using sgpp::base::SPolyBase;
using sgpp::base::SPolyBoundaryBase;

using sgpp::base::SPolyModifiedBase;

using sgpp::base::OperationEval;

BOOST_AUTO_TEST_SUITE(testQuadratureLinear)

BOOST_AUTO_TEST_CASE(testQuadratureLinear) {
  size_t dim = 3;
  size_t level = 4;

  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 5.0));
  grid->getBoundingBox().setBoundary(1, BoundingBox1D(-2.0, 2.0));
  grid->getBoundingBox().setBoundary(2, BoundingBox1D(2.0, 3.0));
  const double determinant = 2.0 * 4.0 * 1.0;  // volume of BoundingBox

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  double qres = 0.0;
  double lSum;

  for (size_t i = 0; i < N; i++) {
    lSum = static_cast<double>(gS.getPoint(i).getLevelSum());
    qres += pow(2, -lSum) * alpha->get(i);
  }

  // take BoundingBox into account by multiplying with determinant
  qres *= determinant;

  double quadOperation =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, qres, 0.0);

  delete alpha;
}

BOOST_AUTO_TEST_CASE(testQuadratureModLinear) {
  size_t dim = 3;
  size_t level = 4;

  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 5.0));
  grid->getBoundingBox().setBoundary(1, BoundingBox1D(-2.0, 2.0));
  grid->getBoundingBox().setBoundary(2, BoundingBox1D(2.0, 3.0));
  const double determinant = 2.0 * 4.0 * 1.0;  // volume of BoundingBox

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  double qres = 0.0;
  double lSum;

  // double resDirect = sgpp::op_factory::createOperationQuadrature(*grid)->doQuadrature(*alpha);

  for (size_t i = 0; i < N; i++) {
    HashGridPoint& gp = gS.getPoint(i);
    lSum = static_cast<double>(gp.getLevelSum());

    // account for the fact that the modified basis functions (left-most and right-most)
    // have integral 2*2^(-lvl) = 2^(-lvl+1)
    for (size_t d = 0; d < dim; d++) {
      HashGridPoint::level_type lvl = gp.getLevel(d);
      HashGridPoint::index_type idx = gp.getIndex(d);
      if ((idx == 1) || (idx == (static_cast<HashGridPoint::index_type>(1) << lvl) - 1)) {
        lSum -= 1.0;
      }
    }

    qres += pow(2, -lSum) * alpha->get(i);
  }

  // take BoundingBox into account by multiplying with determinant
  qres *= determinant;

  double quadOperation =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, qres, 0.0);

  delete alpha;
}

BOOST_AUTO_TEST_CASE(testQuadraturePolyBasis) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 4;

  std::unique_ptr<Grid> grid(Grid::createPolyGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 45.0));
  const double determinant = 42.0;  // volume of BoundingBox

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  SPolyBase* basis = new SPolyBase(deg);

  double quad[7] = {0.666667, 0.333333, 0.333333, 0.168254, 0.164444, 0.164444, 0.168254};
  double quadManual = 0.0;

  for (size_t i = 0; i < N; i++) {
    HashGridPoint& gp = gS.getPoint(i);
    HashGridPoint::level_type lvl = gp.getLevel(0);
    HashGridPoint::index_type idx = gp.getIndex(0);
    double quadSGPP =
        basis->getIntegral(static_cast<unsigned int>(lvl), static_cast<unsigned int>(idx));
    BOOST_CHECK_CLOSE(quadSGPP, quad[i], 5);
    quadManual += alpha->get(i) * quad[i];
  }

  // take BoundingBox into account by multiplying with determinant
  quadManual *= determinant;

  double quadOperation =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, quadManual, 0.00012);

  delete alpha;
  delete basis;
}

BOOST_AUTO_TEST_CASE(testQuadratureModPolyBasis) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 4;

  std::unique_ptr<Grid> grid(Grid::createModPolyGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 45.0));
  const double determinant = 42.0;  // volume of BoundingBox

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  SPolyModifiedBase* basis = new SPolyModifiedBase(deg);

  // correct values
  double quad[7] = {1.0, 0.5, 0.5, 0.25, 0.164444, 0.164444, 0.25};
  double quadManual = 0.0;

  for (size_t i = 0; i < N; i++) {
    HashGridPoint& gp = gS.getPoint(i);
    HashGridPoint::level_type lvl = gp.getLevel(0);
    HashGridPoint::index_type idx = gp.getIndex(0);
    double quadSGPP =
        basis->getIntegral(static_cast<unsigned int>(lvl), static_cast<unsigned int>(idx));
    BOOST_CHECK_CLOSE(quadSGPP, quad[i], 5);
    quadManual += alpha->get(i) * quad[i];
  }

  // take BoundingBox into account by multiplying with determinant
  quadManual *= determinant;

  double quadOperation =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);

  BOOST_CHECK_CLOSE(quadOperation, quadManual, 0.00012);

  delete alpha;
  delete basis;
}

BOOST_AUTO_TEST_CASE(testQuadraturePolyBoundaryBasis) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 4;

  std::unique_ptr<Grid> grid(Grid::createPolyBoundaryGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 45.0));
  const double determinant = 42.0;  // volume of BoundingBox

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  SPolyBoundaryBase* basis = new SPolyBoundaryBase(deg);

  double quad[9] = {0.5, 0.5, 0.666667, 0.333333, 0.333333, 0.168254, 0.164444, 0.164444, 0.168254};
  double quadManual = 0.0;

  for (size_t i = 0; i < N; i++) {
    HashGridPoint& gp = gS.getPoint(i);
    HashGridPoint::level_type lvl = gp.getLevel(0);
    HashGridPoint::index_type idx = gp.getIndex(0);
    double quadSGPP =
        basis->getIntegral(static_cast<unsigned int>(lvl), static_cast<unsigned int>(idx));
    BOOST_CHECK_CLOSE(quadSGPP, quad[i], 5);
    quadManual += alpha->get(i) * quad[i];
  }

  // take BoundingBox into account by multiplying with determinant
  quadManual *= determinant;

  double quadOperation =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);

  BOOST_CHECK_CLOSE(quadOperation, quadManual, 0.0001);

  delete alpha;
  delete basis;
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(testQuadratureMC)

BOOST_AUTO_TEST_CASE(testQuadratureMC) {
  size_t dim = 3;
  size_t level = 3;

  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 5.0));
  grid->getBoundingBox().setBoundary(1, BoundingBox1D(-2.0, 2.0));
  grid->getBoundingBox().setBoundary(2, BoundingBox1D(2.0, 3.0));

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector* alpha = new DataVector(v);

  double resDirect =
      std::unique_ptr<OperationQuadrature>(sgpp::op_factory::createOperationQuadrature(*grid))
          ->doQuadrature(*alpha);

  // Monte Carlo quadrature
  OperationQuadratureMC* opMC = new OperationQuadratureMC(*grid, 100000);
  double resMC = opMC->doQuadrature(*alpha);

  std::cout << resDirect << " " << resMC << std::endl;

  BOOST_CHECK_CLOSE(resDirect, resMC, 1.0);

  delete alpha;
  delete opMC;
}

BOOST_AUTO_TEST_CASE(test_GaussQuadrature) {
  sgpp::base::GaussLegendreQuadRule1D quadRule1D;

  size_t level = 10;
  DataVector pweight(level);
  DataVector coordinates(level);
  quadRule1D.getLevelPointsAndWeights(level, coordinates, pweight);
  std::cout << level << std::endl;
  std::cout << coordinates.toString() << std::endl;
  std::cout << pweight.toString() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(testQuadratureBSpline)

double calc_equidistant_sum(Grid& grid, DataVector& alpha) {
  std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEvalNaive(grid));
  double sum = 0.0;
  double resolution = 10000.0;
  DataVector point(1);
  for (int i = 0; i < static_cast<int>(resolution); i++) {
    point[0] = i / resolution;
    sum += opEval->eval(alpha, point);
  }
  sum /= resolution;
  return sum;
}

BOOST_AUTO_TEST_CASE(testQuadratureBSpline) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  std::unique_ptr<OperationQuadrature> quadOperation(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double quadValue = quadOperation->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadValue, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureModBSpline) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createModBsplineGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  double quadOperation = sgpp::op_factory::createOperationQuadrature(*grid)->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadOperation, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureBSplineBoundary) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineBoundaryGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  std::unique_ptr<OperationQuadrature> quadOperation(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double quadValue = quadOperation->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadValue, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureBSplineClenshawCurtis) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineClenshawCurtisGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  std::unique_ptr<OperationQuadrature> quadOperation(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double quadValue = quadOperation->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadValue, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureModBSplineClenshawCurtis) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineClenshawCurtisGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  std::unique_ptr<OperationQuadrature> quadOperation(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double quadValue = quadOperation->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadValue, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureFundamentalSpline) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineClenshawCurtisGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  std::unique_ptr<OperationQuadrature> quadOperation(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double quadValue = quadOperation->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadValue, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_CASE(testQuadratureModFundamentalSpline) {
  size_t dim = 1;
  size_t level = 3;
  size_t deg = 3;

  std::unique_ptr<Grid> grid(Grid::createBsplineClenshawCurtisGrid(dim, deg));
  grid->getGenerator().regular(level);
  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();
  std::vector<double> v;

  for (size_t i = 0; i < N; i++) v.push_back(static_cast<double>(i));

  DataVector alpha(v);

  double quadOperation = sgpp::op_factory::createOperationQuadrature(*grid)->doQuadrature(alpha);

  double equidistant_sum = calc_equidistant_sum(*grid, alpha);

  BOOST_CHECK_CLOSE(quadOperation, equidistant_sum, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
