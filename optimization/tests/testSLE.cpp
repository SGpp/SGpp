#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/InterpolantFunction.hpp>
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
      A.set(i, j, Aij);
      Ax[i] += Aij * x.get(j);

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
    SGPP::float_t ri = b.get(i);

    for (size_t j = 0; j < n; j++) {
      ri -= A.get(i, j) * x.get(j);
    }

    rNormSquared += ri * ri;
    bNormSquared += b.get(i) * b.get(i);
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

  // test various SLE dimensions
  for (size_t n : {
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200
       }) {
    // generate random matrix and RHS
    base::DataMatrix A(n, n);
    base::DataVector b(n);

    for (size_t i = 0; i < n; i++) {
      b[i] = randomNumberGenerator.getUniformRN(-1.0, 1.0);

      for (size_t j = 0; j < n; j++) {
        A.set(i, j, randomNumberGenerator.getUniformRN(-1.0, 1.0));
      }
    }

    // full SLE
    FullSLE system(A);

    for (const auto& solver : solvers) {
      if ((dynamic_cast<sle_solver::BiCGStab*>(solver.get()) != nullptr) &&
          (n > (use_double_precision ? 20 : 8))) {
        /*
         * BiCGStab is really weak and can't solve bigger systems
         * (a bug in the implementation is unlikely as MATLAB
         * shows the same result, but BiCGStab should only be used for Newton's
         * method of optimization - for hierarchisation, only external solvers
         * should be used)
         */
        continue;
      } else if ((dynamic_cast<sle_solver::Gmmpp*>(solver.get()) != nullptr) &&
                 (!use_double_precision) && (n == 200)) {
        // Gmm++ doesn't converge using single precision for larger systems
        continue;
      }

      base::DataVector x(0);

      // solve system
      BOOST_CHECK(solver->solve(system, b, x));

      // test solution
      testSLESolution(A, x, b);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestHierarchization) {
  // Test SGPP::optimization::HierarchisationSLE.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t l = 4;

  base::DataVector x(d);
  ExampleFunction f;
  sle_solver::Auto solver;

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
    // generate regular sparse grid
    std::unique_ptr<base::GridGenerator> gridGen(grid->createGridGenerator());
    gridGen->regular(l);
    const size_t n = grid->getSize();
    base::DataVector functionValues(n);

    for (size_t i = 0; i < n; i++) {
      base::GridIndex& gp = *grid->getStorage()->get(i);

      // don't forget to set the point distribution to Clenshaw-Curtis
      // if necessary (currently not done automatically)
      if ((std::string(grid->getType()) == "bsplineClenshawCurtis") ||
          (std::string(grid->getType()) == "linearClenshawCurtis")) {
        gp.setPointDistribution(
          base::GridIndex::PointDistribution::ClenshawCurtis);
      }

      for (size_t t = 0; t < d; t++) {
        x[t] = gp.getCoord(t);
      }

      functionValues[i] = f.eval(x);
    }

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
    InterpolantFunction ft(*grid, alpha);

    for (size_t i = 0; i < 100; i++) {
      for (size_t t = 0; t < d; t++) {
        // don't go near the boundary (should suffice)
        x[t] = randomNumberGenerator.getUniformRN(0.2, 0.8);
      }

      // test infinity norm of difference roughly
      BOOST_CHECK_SMALL(f.eval(x) - ft.eval(x),
                        static_cast<SGPP::float_t>(0.3));
    }
  }
}
