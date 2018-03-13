// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/globaldef.hpp>
namespace sgpp {
namespace pde {

double uniform_distributed_approximation(sgpp::base::Grid& grid, size_t i, size_t j) {
  const size_t d = grid.getDimension();
  const size_t resolution = 10000;
  const double h = 1.0 / static_cast<double>(resolution);
  sgpp::base::GridStorage& storage = grid.getStorage();
  sgpp::base::SBasis& basis = grid.getBasis();
  double res = 1.0;
  for (size_t k = 0; k < d; k++) {
    const sgpp::base::level_t lik = storage.getPoint(i).getLevel(k);
    const sgpp::base::level_t ljk = storage.getPoint(j).getLevel(k);
    const sgpp::base::index_t iik = storage.getPoint(i).getIndex(k);
    const sgpp::base::index_t ijk = storage.getPoint(j).getIndex(k);
    // trapezoidal rule
    double temp_res = 0.0;
    // --------------------------------------------------------------------------
    // apply trapezoidal rule
    temp_res += basis.eval(lik, iik, 0.0) * basis.eval(ljk, ijk, 0.0) / 2.;
    for (size_t c = 1; c < resolution; c++) {
      double x = static_cast<double>(c) * h;
      temp_res += basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
    }
    temp_res += basis.eval(lik, iik, 1.0) * basis.eval(ljk, ijk, 1.0) / 2.;
    // --------------------------------------------------------------------------
    res *= temp_res * h;
  }
  return res;
}

BOOST_AUTO_TEST_SUITE(testOperationMatrixLTwoDotExplicit)

// test for Linear
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitLinear) {
  const size_t d = 3;
  const size_t l = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createLinearGrid(d));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for ModLinear
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitModLinear) {
  const size_t d = 3;
  const size_t l = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createModLinearGrid(d));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for Poly
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitPoly) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createPolyGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for PolyBoundary
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitPolyBoundary) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createPolyBoundaryGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for ModPoly
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitModPoly) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createModPolyGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for PolyClenshawCurtis
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitPolyClenshawCurtis) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for PolyClenshawCurtisBoundary
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary) {
  const size_t d = 2;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createPolyClenshawCurtisBoundaryGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if  explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for ModPolyClenshawCurtis
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitModPolyClenshawCurtis) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createModPolyClenshawCurtisGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  // test Explicit correctness
  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for Bspline
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitBspline) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for BsplineBoundary
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitBsplineBoundary) {
  const size_t d = 2;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineBoundaryGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for ModBspline
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitModBspline) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createModBsplineGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

// test for BsplineClenshawCurtis
BOOST_AUTO_TEST_CASE(testOperationMatrixLTwoDotExplicitBsplineClenshawCurtis) {
  const size_t d = 3;
  const size_t l = 3;
  const size_t p = 3;
  sgpp::base::Grid* grid(sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, p));
  grid->getGenerator().regular(l);

  sgpp::base::DataMatrix m(grid->getSize(), grid->getSize());
  sgpp::base::OperationMatrix* opExplicit =
      sgpp::op_factory::createOperationLTwoDotExplicit(&m, *grid);

  for (size_t i = 0; i < grid->getSize(); i++) {
    for (size_t j = i; j < grid->getSize(); j++) {
      double approx = uniform_distributed_approximation(*grid, i, j);
      // std::cout << std::abs(approx - m.get(i, j)) << std::endl;
      BOOST_CHECK_SMALL(approx - m.get(i, j), 1e-3);
    }
  }

  // test if explicit equals implicit
  sgpp::base::OperationMatrix* opImplicit = sgpp::op_factory::createOperationLTwoDotProduct(*grid);
  sgpp::base::DataVector alpha(grid->getSize(), 1.0);

  sgpp::base::DataVector resultImplicit(grid->getSize());
  sgpp::base::DataVector resultExplicit(grid->getSize());

  opExplicit->mult(alpha, resultExplicit);
  opImplicit->mult(alpha, resultImplicit);
  for (size_t i = 0; i < grid->getSize(); i++) {
    BOOST_CHECK_SMALL(resultImplicit.get(i) - resultExplicit.get(i), 1e-12);
  }

  delete grid;
  delete opImplicit;
  delete opExplicit;
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace pde
}  // namespace sgpp
