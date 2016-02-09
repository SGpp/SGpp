#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>


using namespace SGPP::base;

BOOST_AUTO_TEST_SUITE(testQuadratureLinear)

BOOST_AUTO_TEST_CASE(testQuadrature) {

  size_t dim = 3;
  size_t level = 4;

  Grid* grid = Grid::createLinearGrid(dim);

  grid->createGridGenerator()->regular(level);
  GridStorage* gS = grid->getStorage();

  size_t N = gS->size();
  std::vector<int> v(N);

  for (size_t i = 0; i < N; i++)
    v.push_back(static_cast<int>(i));

  DataVector* alpha = new DataVector(v);

  SGPP::float_t qres = 0.0;
  SGPP::float_t lSum;

  for (size_t i = 0; i < N; i++) {
    lSum = static_cast<SGPP::float_t>(gS->get(i)->getLevelSum());
    qres += pow(2, -lSum) * alpha->get(i);
  }

  OperationQuadrature* quadOp = SGPP::op_factory::createOperationQuadrature(
                                  *grid);
  SGPP::float_t quadOperation = quadOp->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, qres, 0.0);

  delete grid;
  delete alpha;
  delete quadOp;
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(testQuadratureMC)

BOOST_AUTO_TEST_CASE(testQuadratureMC) {

  size_t dim = 3;
  size_t level = 3;

  Grid* grid = Grid::createLinearGrid(dim);

  grid->createGridGenerator()->regular(level);
  GridStorage* gS = grid->getStorage();
  size_t N = gS->size();
  std::vector<int> v(N);

  for (size_t i = 0; i < N; i++)
    v[static_cast<int>(i)] = static_cast<int>(i);

  DataVector* alpha = new DataVector(v);

  OperationQuadrature* quadOp = SGPP::op_factory::createOperationQuadrature(
                                  *grid);
  SGPP::float_t resDirect = quadOp->doQuadrature(*alpha);

  // Monte Carlo quadrature
  OperationQuadratureMC* opMC = new OperationQuadratureMC(*grid, 100000);
  SGPP::float_t resMC = opMC->doQuadrature(*alpha);

  BOOST_CHECK_CLOSE(resDirect, resMC, 1.0);

  delete grid;
  delete alpha;
  delete opMC;
  delete quadOp;
}

BOOST_AUTO_TEST_CASE(testQuadraturePolyBasis) {

  size_t dim = 1;
  size_t level = 3;
  size_t deg = 4;

  Grid* grid = Grid::createPolyGrid(dim, deg);
  grid->createGridGenerator()->regular(level);
  GridStorage* gS = grid->getStorage();
  size_t N = gS->size();
  std::vector<int> v(N);

  for (size_t i = 0; i < N; i++)
    v[static_cast<int>(i)] = static_cast<int>(i);

  DataVector* alpha = new DataVector(v);

  SPolyBase* basis = new SPolyBase(deg);

  SGPP::float_t quad[7] = {0.666667, 0.333333, 0.333333, 0.168254, 0.164444, 0.164444, 0.168254};
  SGPP::float_t quadManual = 0.0;

  for (size_t i = 0; i < N; i++) {
    HashGridIndex* gp = gS->get(i);
    HashGridIndex::level_type lvl = gp->getLevel(0);
    HashGridIndex::index_type idx = gp->getIndex(0);
    SGPP::float_t quadSGPP = basis->getIntegral(static_cast<unsigned int>(lvl),
                             static_cast<unsigned int>(idx));
    BOOST_CHECK_CLOSE(quadSGPP, quad[i], 5);
    quadManual += alpha->get(i) * quad[i];
  }

  OperationQuadrature* quadOp = SGPP::op_factory::createOperationQuadrature(
                                  *grid);
  SGPP::float_t quadOperation = quadOp->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, quadManual, 0.00012);

  delete grid;
  delete alpha;
  delete basis;
  delete quadOp;
}
BOOST_AUTO_TEST_CASE(testQuadraturePolyBoundaryBasis) {

  size_t dim = 1;
  size_t level = 3;
  size_t deg = 4;

  Grid* grid = Grid::createPolyBoundaryGrid(dim, deg);
  grid->createGridGenerator()->regular(level);
  GridStorage* gS = grid->getStorage();
  size_t N = gS->size();
  std::vector<int> v(N);

  for (size_t i = 0; i < N; i++)
    v[static_cast<int>(i)] = static_cast<int>(i);

  DataVector* alpha = new DataVector(v);

  SPolyBoundaryBase* basis = new SPolyBoundaryBase(deg);

  SGPP::float_t quad[9] = {0.5, 0.5, 0.666667, 0.333333, 0.333333, 0.168254, 0.164444, 0.164444, 0.168254};
  SGPP::float_t quadManual = 0.0;

  for (size_t i = 0; i < N; i++) {
    HashGridIndex* gp = gS->get(i);
    HashGridIndex::level_type lvl = gp->getLevel(0);
    HashGridIndex::index_type idx = gp->getIndex(0);
    SGPP::float_t quadSGPP = basis->getIntegral(static_cast<unsigned int>(lvl),
                             static_cast<unsigned int>(idx));
    BOOST_CHECK_CLOSE(quadSGPP, quad[i], 5);
    quadManual += alpha->get(i) * quad[i];
  }

  OperationQuadrature* quadOp = SGPP::op_factory::createOperationQuadrature(
                                  *grid);
  SGPP::float_t quadOperation = quadOp->doQuadrature(*alpha);
  BOOST_CHECK_CLOSE(quadOperation, quadManual, 0.0001);

  delete grid;
  delete alpha;
  delete basis;
  delete quadOp;
}
BOOST_AUTO_TEST_SUITE_END()
