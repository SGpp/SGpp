// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>

class ExampleFunction : public sgpp::optimization::ScalarFunction {
 public:
  ExampleFunction() : sgpp::optimization::ScalarFunction(2) {}

  double eval(const sgpp::base::DataVector& x) {
    // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
    return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  }

  virtual void clone(std::unique_ptr<sgpp::optimization::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunction>(new ExampleFunction(*this));
  }
};

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  sgpp::optimization::Printer::getInstance().setVerbosity(-1);

  // objective function
  ExampleFunction f;
  // dimension of domain
  const size_t d = f.getNumberOfParameters();
  // B-spline degree
  const size_t p = 3;
  // maximal number of grid points
  const size_t N = 200;
  // adaptivity of grid generation
  const double gamma = 0.95;

  /**
   * First, we define a grid with modified B-spline basis functions and
   * an iterative grid generator, which can generate the grid adaptively.
   */
  sgpp::base::ModBsplineGrid grid(d, p);
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(f, grid, N, gamma);

  /**
   * With the iterative grid generator, we generate adaptively a sparse grid.
   */

  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  /**
   * Then, we hierarchize the function values to get hierarchical B-spline
   * coefficients of the B-spline sparse grid interpolant
   * \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
   */
  sgpp::base::DataVector functionValues(gridGen.getFunctionValues());
  sgpp::base::DataVector coeffs(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(grid);
  sgpp::optimization::sle_solver::UMFPACK UMFPACKSolver;

  sgpp::combigrid::Stopwatch watch;
  watch.start();
  if (!UMFPACKSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }
  std::cout << "UMFPACK: " << watch.elapsedSeconds() << std::endl;

  sgpp::optimization::sle_solver::Armadillo ArmadilloSolver;

  watch.start();
  if (!ArmadilloSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }
  std::cout << "Armadillo: " << watch.elapsedSeconds() << std::endl;

  sgpp::optimization::sle_solver::Auto AutoSolver;

  watch.start();
  if (!AutoSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }
  std::cout << "Auto: " << watch.elapsedSeconds() << std::endl;

  return 0;
}
