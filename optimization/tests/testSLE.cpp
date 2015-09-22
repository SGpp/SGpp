#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/sle/solver/Gmmpp.hpp>
#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include "ObjectiveFunctions.hpp"
#include "GridCreator.hpp"

const bool use_double_precision =
#if USE_DOUBLE_PRECISION
  true;
#else
  false;
#endif /* USE_DOUBLE_PRECISION */

using namespace SGPP;
using namespace SGPP::optimization;

void testSLESystem(SLE& system, const base::DataVector& x,
                   const base::DataVector& b,
                   base::DataMatrix& A) {
  // Test SGPP::optimization::SLE::getMatrixEntry, isMatrixEntryNonZero and
  // matrixVectorMultiplication. Returns system matrix as pysgpp.DataMatrix.
  const size_t n = x.getSize();
  BOOST_CHECK_EQUAL(system.getDimension(), n);
  A.resize(n, n);
  base::DataVector Ax(n, 0.0);

  // A*x calculated directly
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      const SGPP::float_t Aij = system.getMatrixEntry(i, j);
      A(i, j) = Aij;
      Ax[i] += Aij * x[j];

      // test isMatrixEntryNonZero
      BOOST_CHECK_EQUAL(system.isMatrixEntryNonZero(i, j), Aij != 0);
    }
  }

  // A*x calculated by SGPP::optimization
  base::DataVector Ax2(0);
  system.matrixVectorMultiplication(x, Ax2);

  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_CLOSE(Ax[i], Ax2[i], 1e-10);
  }
}

void testSLESolution(const base::DataMatrix& A,
                     base::DataVector& x,
                     const base::DataVector& b) {
  const size_t n = b.getSize();
  BOOST_CHECK_EQUAL(x.getSize(), n);
  SGPP::float_t rNormSquared = 0.0;
  SGPP::float_t bNormSquared = 0.0;

  for (size_t i = 0; i < n; i++) {
    SGPP::float_t ri = b[i];
    bNormSquared += ri * ri;

    for (size_t j = 0; j < n; j++) {
      ri -= A(i, j) * x[j];
    }

    rNormSquared += ri * ri;
  }

  // test relative residual
  BOOST_CHECK_SMALL(std::sqrt(rNormSquared / bNormSquared),
                    (use_double_precision ? static_cast<SGPP::float_t>(1e-6) :
                     static_cast<SGPP::float_t>(1e-3)));
}

BOOST_AUTO_TEST_CASE(TestSLESolvers) {
  // Test SGPP::optimization::sle_solver with SGPP::optimization::FullSLE.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t m = 4;

  // default solvers
  std::vector<std::unique_ptr<sle_solver::SLESolver>> solvers;
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::BiCGStab())));
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::GaussianElimination())));
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::Auto())));

  // additional solvers if SGPP::opt was compiled with them
#ifdef USE_ARMADILLO
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::Armadillo())));
#endif /* USE_ARMADILLO */
#ifdef USE_EIGEN
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::Eigen())));
#endif /* USE_EIGEN */
#ifdef USE_GMMPP
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::Gmmpp())));
#endif /* USE_GMMPP */
#ifdef USE_UMFPACK
  solvers.push_back(std::move(std::unique_ptr<sle_solver::SLESolver>(
                                new sle_solver::UMFPACK())));
#endif /* USE_UMFPACK */

  // test getters/setters
  {
    sle_solver::BiCGStab biCGStab;

    const size_t maxItCount = 42;
    biCGStab.setMaxItCount(maxItCount);
    BOOST_CHECK_EQUAL(biCGStab.getMaxItCount(), maxItCount);

    const SGPP::float_t tolerance = 0.42;
    biCGStab.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(biCGStab.getTolerance(), tolerance);

    base::DataVector startingPoint(3);
    startingPoint[0] = 1.2;
    startingPoint[1] = 3.4;
    startingPoint[2] = 5.6;
    biCGStab.setStartingPoint(startingPoint);
    base::DataVector startingPoint2 = biCGStab.getStartingPoint();
    BOOST_CHECK_EQUAL(startingPoint.getSize(), startingPoint2.getSize());

    for (size_t t = 0; t < startingPoint.getSize(); t++) {
      BOOST_CHECK_EQUAL(startingPoint[t], startingPoint2[t]);
    }
  }

  // test various SLE dimensions
  for (size_t n : {
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200
       }) {
    // generate random matrix and RHS
    base::DataMatrix A(n, n);
    base::DataVector b(n);
    base::DataMatrix B(n, m);

    for (size_t i = 0; i < n; i++) {
      b[i] = randomNumberGenerator.getUniformRN(-1.0, 1.0);

      for (size_t j = 0; j < m; j++) {
        B(i, j) = randomNumberGenerator.getUniformRN(-1.0, 1.0);
      }

      for (size_t j = 0; j < n; j++) {
        A(i, j) = randomNumberGenerator.getUniformRN(-1.0, 1.0);
      }
    }

    // full SLE
    FullSLE system(A);

    for (const auto& solver : solvers) {
      if ((dynamic_cast<sle_solver::BiCGStab*>(solver.get()) != nullptr) &&
          (n > (use_double_precision ? 20 : 6))) {
        /*
         * BiCGStab is really weak and can't solve bigger systems
         * (a bug in the implementation is unlikely as MATLAB
         * shows the same result, but BiCGStab should only be used for Newton's
         * method of optimization - for hierarchisation, only external solvers
         * should be used)
         */
        continue;
      } else if ((dynamic_cast<sle_solver::Gmmpp*>(solver.get()) != nullptr) &&
                 (!use_double_precision) && (n > 50)) {
        // Gmm++ doesn't converge using single precision for larger systems
        continue;
      }

      // solve system and test solution
      base::DataVector x(0);
      BOOST_CHECK(solver->solve(system, b, x));
      testSLESolution(A, x, b);

      base::DataMatrix X(0, 0);
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
  // Test SGPP::optimization::FullSLE.
  base::DataMatrix A(3, 3, 0.0);
  A(0, 1) = 12.3;
  A(1, 2) = 42.1337;

  FullSLE sle(A);
  std::unique_ptr<CloneableSLE> sle2(nullptr);
  sle.clone(sle2);

  BOOST_CHECK_EQUAL(sle2->countNNZ(), 2);

  if (use_double_precision) {
    BOOST_CHECK_EQUAL(sle2->getMatrixEntry(0, 1), 12.3);
    BOOST_CHECK_EQUAL(sle2->getMatrixEntry(1, 2), 42.1337);
  } else {
    BOOST_CHECK_CLOSE(sle2->getMatrixEntry(0, 1), 12.3, 1e-5);
    BOOST_CHECK_CLOSE(sle2->getMatrixEntry(1, 2), 42.1337, 1e-5);
  }

  BOOST_CHECK(sle2->isMatrixEntryNonZero(1, 2));
  BOOST_CHECK(!sle2->isMatrixEntryNonZero(2, 2));
}

BOOST_AUTO_TEST_CASE(TestHierarchisationSLE) {
  // Test SGPP::optimization::HierarchisationSLE.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t l = 4;

  base::DataVector x(d);
  sle_solver::Auto solver;
  ExampleFunction f;

  // Test All The Grids!
  std::vector<std::unique_ptr<base::Grid>> grids;
  createSupportedGrids(d, p, grids);

  for (auto& grid : grids) {
    base::DataVector functionValues(0);
    createSampleGrid(*grid, l, f, functionValues);

    // create hierarchization system
    HierarchisationSLE system(*grid);
    base::DataVector alpha(0);

    // solve system
    BOOST_CHECK(solver.solve(system, functionValues, alpha));

    // test system
    base::DataMatrix A(0, 0);
    testSLESystem(system, alpha, functionValues, A);

    // test solution
    testSLESolution(A, alpha, functionValues);

    // create interpolant
    InterpolantScalarFunction ft(*grid, alpha);

    // test InterpolantScalarFunction::clone()
    std::unique_ptr<ScalarFunction> ft2(nullptr);
    ft.clone(ft2);

    for (size_t i = 0; i < 100; i++) {
      for (size_t t = 0; t < d; t++) {
        // don't go near the boundary (should suffice)
        x[t] = randomNumberGenerator.getUniformRN(0.2, 0.8);
      }

      // test infinity norm of difference roughly
      BOOST_CHECK_SMALL(f.eval(x) - ft2->eval(x),
                        static_cast<SGPP::float_t>(0.3));
    }
  }
}
