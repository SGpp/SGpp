// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/solver/BiCGStab.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/base/tools/sle/solver/Gmmpp.hpp>
#include <sgpp/base/tools/sle/solver/UMFPACK.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>

#include <cmath>
#include <vector>

using sgpp::base::CloneableSLE;
using sgpp::base::FullSLE;
using sgpp::base::HierarchisationSLE;
using sgpp::base::InterpolantScalarFunction;
using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;
using sgpp::base::ScalarFunction;
using sgpp::base::SLE;

// The following example function and createGrid routines are according to the routines in
// sgpp::optimization GridCreator and ObjectiveFunctions. Because the SLE routines were moved to
// base, these functionalities were rewritten here to be independent of optimization

class ExampleFunctionSLE : public ScalarFunction {
 public:
  ExampleFunctionSLE() : ScalarFunction(2) {}
  ~ExampleFunctionSLE() override {}

  double eval(const sgpp::base::DataVector& x) override {
    // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
    if ((x[0] >= 0.0) && (x[0] <= 1.0) && (x[1] >= 0.0) && (x[1] <= 1.0)) {
      return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
    } else {
      return INFINITY;
    }
  }
  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new ExampleFunctionSLE(*this));
  }
};

void createSupportedGridsSLE(size_t d, size_t p,
                             std::vector<std::unique_ptr<sgpp::base::Grid>>& grids) {
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createBsplineGrid(d, p)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createBsplineBoundaryGrid(d, p)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createModBsplineClenshawCurtisGrid(d, p)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createFundamentalSplineGrid(d, p)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModFundamentalSplineGrid(d, p)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(
      sgpp::base::Grid::createLinearClenshawCurtisBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createWaveletGrid(d)));
  grids.push_back(
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createWaveletBoundaryGrid(d)));
  grids.push_back(std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModWaveletGrid(d)));
}

void createSampleGridSLE(sgpp::base::Grid& grid, size_t l, sgpp::base::ScalarFunction& f,
                         sgpp::base::DataVector& functionValues) {
  sgpp::base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  // generate regular sparse grid
  gridStorage.clear();
  grid.getGenerator().regular(l);
  const size_t n = gridStorage.getSize();
  sgpp::base::DataVector x(d);

  functionValues.resize(n);

  for (size_t i = 0; i < n; i++) {
    x = gridStorage.getCoordinates(gridStorage[i]);
    functionValues[i] = f.eval(x);
  }
}

void testSLESystem(SLE& system, const sgpp::base::DataVector& x, const sgpp::base::DataVector& b,
                   sgpp::base::DataMatrix& A) {
  // Test sgpp::base::SLE::getMatrixEntry, isMatrixEntryNonZero and
  // matrixVectorMultiplication. Returns system matrix as pysgpp.DataMatrix.
  const size_t n = x.getSize();
  BOOST_CHECK_EQUAL(system.getDimension(), n);
  A.resize(n, n);
  sgpp::base::DataVector Ax(n, 0.0);

  // A*x calculated directly
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      const double Aij = system.getMatrixEntry(i, j);
      A(i, j) = Aij;
      Ax[i] += Aij * x[j];

      // test isMatrixEntryNonZero
      BOOST_CHECK_EQUAL(system.isMatrixEntryNonZero(i, j), Aij != 0);
    }
  }

  // A*x calculated by sgpp::optimization
  sgpp::base::DataVector Ax2(0);
  system.matrixVectorMultiplication(x, Ax2);

  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_CLOSE(Ax[i], Ax2[i], 1e-10);
  }
}

void testSLESolution(const sgpp::base::DataMatrix& A, sgpp::base::DataVector& x,
                     const sgpp::base::DataVector& b) {
  const size_t n = b.getSize();
  BOOST_CHECK_EQUAL(x.getSize(), n);
  double rNormSquared = 0.0;
  double bNormSquared = 0.0;

  for (size_t i = 0; i < n; i++) {
    double ri = b[i];
    bNormSquared += ri * ri;

    for (size_t j = 0; j < n; j++) {
      ri -= A(i, j) * x[j];
    }

    rNormSquared += ri * ri;
  }

  // test relative residual
  BOOST_CHECK_SMALL(std::sqrt(rNormSquared / bNormSquared), 1e-6);
}

BOOST_AUTO_TEST_CASE(TestSLESolvers) {
  // Test sgpp::base::sle_solver with sgpp::base::FullSLE.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t m = 4;

  // default solvers
  std::vector<std::unique_ptr<sgpp::base::sle_solver::SLESolver>> solvers;
  solvers.push_back(
      std::unique_ptr<sgpp::base::sle_solver::SLESolver>(new sgpp::base::sle_solver::BiCGStab()));
  solvers.push_back(std::unique_ptr<sgpp::base::sle_solver::SLESolver>(
      new sgpp::base::sle_solver::GaussianElimination()));
  solvers.push_back(
      std::unique_ptr<sgpp::base::sle_solver::SLESolver>(new sgpp::base::sle_solver::Auto()));

  // additional solvers if sgpp::opt was compiled with them
#ifdef USE_ARMADILLO
  solvers.push_back(
      std::unique_ptr<sgpp::base::sle_solver::SLESolver>(new sgpp::base::sle_solver::Armadillo()));
#endif /* USE_ARMADILLO */
#ifdef USE_EIGEN
  solvers.push_back(
      std::unique_ptr<sgpp::base::sle_solver::SLESolver>(new sgpp::base::sle_solver::Eigen()));
#endif /* USE_EIGEN */
#ifdef USE_GMMPP
  solvers.push_back(std::unique_ptr<sgpp::base::sle_solver::SLESolver>(
      new sgpp::base::sle_solver::Gmmpp()));
#endif /* USE_GMMPP */
#ifdef USE_UMFPACK
  solvers.push_back(std::unique_ptr<sgpp::base::sle_solver::SLESolver>(
      new sgpp::base::sle_solver::UMFPACK()));
#endif /* USE_UMFPACK */

  // test getters/setters
  {
    sgpp::base::sle_solver::BiCGStab biCGStab;

    const size_t maxItCount = 42;
    biCGStab.setMaxItCount(maxItCount);
    BOOST_CHECK_EQUAL(biCGStab.getMaxItCount(), maxItCount);

    const double tolerance = 0.42;
    biCGStab.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(biCGStab.getTolerance(), tolerance);

    sgpp::base::DataVector startingPoint(3);
    startingPoint[0] = 1.2;
    startingPoint[1] = 3.4;
    startingPoint[2] = 5.6;
    biCGStab.setStartingPoint(startingPoint);
    sgpp::base::DataVector startingPoint2 = biCGStab.getStartingPoint();
    BOOST_CHECK_EQUAL(startingPoint.getSize(), startingPoint2.getSize());

    for (size_t t = 0; t < startingPoint.getSize(); t++) {
      BOOST_CHECK_EQUAL(startingPoint[t], startingPoint2[t]);
    }
  }

  // test various SLE dimensions
  for (size_t n : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200}) {
    // generate random matrix and RHS
    sgpp::base::DataMatrix A(n, n);
    sgpp::base::DataVector b(n);
    sgpp::base::DataMatrix B(n, m);

    for (size_t i = 0; i < n; i++) {
      b[i] = RandomNumberGenerator::getInstance().getUniformRN(-1.0, 1.0);

      for (size_t j = 0; j < m; j++) {
        B(i, j) = RandomNumberGenerator::getInstance().getUniformRN(-1.0, 1.0);
      }

      for (size_t j = 0; j < n; j++) {
        A(i, j) = RandomNumberGenerator::getInstance().getUniformRN(-1.0, 1.0);
      }
    }

    // full SLE
    FullSLE system(A);

    for (const auto& solver : solvers) {
      if ((dynamic_cast<sgpp::base::sle_solver::BiCGStab*>(solver.get()) != nullptr) && (n > 20)) {
        /*
         * BiCGStab is really weak and can't solve bigger systems
         * (a bug in the implementation is unlikely as MATLAB
         * shows the same result, but BiCGStab should only be used for Newton's
         * method of optimization - for hierarchisation, only external solvers
         * should be used)
         */
        continue;
      }

      // solve system and test solution
      sgpp::base::DataVector x(0);
      BOOST_CHECK(solver->solve(system, b, x));
      testSLESolution(A, x, b);

      sgpp::base::DataMatrix X(0, 0);
      BOOST_CHECK(solver->solve(system, B, X));

      for (size_t j = 0; j < m; j++) {
        X.getColumn(j, x);
        B.getColumn(j, b);
        testSLESolution(A, x, b);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestFullSLE) {
  // Test sgpp::base::FullSLE.
  sgpp::base::DataMatrix A(3, 3, 0.0);
  A(0, 1) = 12.3;
  A(1, 2) = 42.1337;

  FullSLE sle(A);
  std::unique_ptr<CloneableSLE> sle2(nullptr);
  sle.clone(sle2);

  BOOST_CHECK_EQUAL(sle2->countNNZ(), 2U);

  BOOST_CHECK_EQUAL(sle2->getMatrixEntry(0, 1), 12.3);
  BOOST_CHECK_EQUAL(sle2->getMatrixEntry(1, 2), 42.1337);

  BOOST_CHECK(sle2->isMatrixEntryNonZero(1, 2));
  BOOST_CHECK(!sle2->isMatrixEntryNonZero(2, 2));
}

BOOST_AUTO_TEST_CASE(TestHierarchisationSLE) {
  // Test sgpp::base::HierarchisationSLE.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t l = 4;

  sgpp::base::DataVector x(d);
  sgpp::base::sle_solver::Auto solver;
  ExampleFunctionSLE f;

  // Test All The Grids!
  std::vector<std::unique_ptr<sgpp::base::Grid>> grids;
  createSupportedGridsSLE(d, p, grids);

  for (auto& grid : grids) {
    sgpp::base::DataVector functionValues(0);
    createSampleGridSLE(*grid, l, f, functionValues);

    // create hierarchization system
    HierarchisationSLE system(*grid);
    sgpp::base::DataVector alpha(0);

    // solve system
    BOOST_CHECK(solver.solve(system, functionValues, alpha));

    // test system
    sgpp::base::DataMatrix A(0, 0);
    testSLESystem(system, alpha, functionValues, A);

    // test solution
    testSLESolution(A, alpha, functionValues);

    // create interpolant
    sgpp::base::InterpolantScalarFunction ft(*grid, alpha);

    // test InterpolantScalarFunction::clone()
    std::unique_ptr<ScalarFunction> ft2(nullptr);
    ft.clone(ft2);

    for (size_t i = 0; i < 100; i++) {
      for (size_t t = 0; t < d; t++) {
        // don't go near the boundary (should suffice)
        x[t] = RandomNumberGenerator::getInstance().getUniformRN(0.2, 0.8);
      }

      // test infinity norm of difference roughly
      BOOST_CHECK_SMALL(f.eval(x) - ft2->eval(x), 0.3);
    }
  }
}
