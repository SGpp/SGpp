// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>

/**
 * Example test function.
 */
class ExampleFunction : public sgpp::optimization::ScalarFunction {
 public:
  /**
   * Constructor.
   */
  ExampleFunction() : sgpp::optimization::ScalarFunction(2) {
  }

  /**
   * Evaluates the test function.
   *
   * @param x     point \f$\vec{x} \in [0, 1]^2\f$
   * @return      \f$f(\vec{x})\f$
   */
  double eval(const sgpp::base::DataVector& x) {
    // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
    return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(
    std::unique_ptr<sgpp::optimization::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunction>(
              new ExampleFunction(*this));
  }
};

/**
 * Prints a separator line.
 */
void printLine() {
  std::cout << "----------------------------------------"
            "----------------------------------------\n";
}

/**
 * Main method.
 *
 * @param argc ignored
 * @param argv ignored
 */
int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "sgpp::optimization example program started.\n\n";
  sgpp::optimization::Printer::getInstance().setVerbosity(2);

  // objective function
  ExampleFunction f;
  // dimension of domain
  const size_t d = f.getNumberOfParameters();
  // B-spline degree
  const size_t p = 3;
  // maximal number of grid points
  const size_t N = 30;
  // adaptivity of grid generation
  const double gamma = 0.95;

  sgpp::base::ModBsplineGrid grid(d, p);
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(
    f, grid, N, gamma);

  // //////////////////////////////////////////////////////////////////////////
  // GRID GENERATION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Generating grid...\n\n";

  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  // //////////////////////////////////////////////////////////////////////////
  // HIERARCHIZATION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Hierarchizing...\n\n";
  sgpp::base::DataVector functionValues(gridGen.getFunctionValues());
  sgpp::base::DataVector coeffs(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(grid);
  sgpp::optimization::sle_solver::Auto sleSolver;

  // solve linear system
  if (!sleSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  // //////////////////////////////////////////////////////////////////////////
  // OPTIMIZATION OF THE SMOOTH INTERPOLANT
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing smooth interpolant...\n\n";
  sgpp::optimization::InterpolantScalarFunction ft(grid, coeffs);
  sgpp::optimization::InterpolantScalarFunctionGradient ftGradient(grid, coeffs);
  sgpp::optimization::optimizer::GradientDescent gradientMethod(ft, ftGradient);
  sgpp::base::DataVector x0(d);
  double fX0;
  double ftX0;

  // determine best grid point as starting point
  {
    sgpp::base::GridStorage& gridStorage = grid.getStorage();

    // index of grid point with minimal function value
    size_t x0Index = std::distance(
                       functionValues.getPointer(),
                       std::min_element(functionValues.getPointer(),
                                        functionValues.getPointer() +
                                        functionValues.getSize()));

    for (size_t t = 0; t < d; t++) {
      x0[t] = gridStorage[x0Index]->getCoord(t);
    }

    fX0 = functionValues[x0Index];
    ftX0 = ft.eval(x0);
  }

  std::cout << "x0 = " << x0.toString() << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  gradientMethod.setStartingPoint(x0);
  gradientMethod.optimize();
  const sgpp::base::DataVector& xOpt = gradientMethod.getOptimalPoint();
  const double ftXOpt = gradientMethod.getOptimalValue();
  const double fXOpt = f.eval(xOpt);

  std::cout << "\nxOpt = " << xOpt.toString() << "\n";
  std::cout << "f(xOpt) = " << fXOpt << ", ft(xOpt) = " << ftXOpt << "\n\n";

  // //////////////////////////////////////////////////////////////////////////
  // NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing objective function (for comparison)...\n\n";

  sgpp::optimization::optimizer::NelderMead nelderMead(f, 1000);
  nelderMead.optimize();
  sgpp::base::DataVector xOptNM = nelderMead.getOptimalPoint();
  const double fXOptNM = nelderMead.getOptimalValue();
  const double ftXOptNM = ft.eval(xOptNM);

  std::cout << "\nnxOptNM = " << xOptNM.toString() << "\n";
  std::cout << "f(xOptNM) = " << fXOptNM <<
            ", ft(xOptNM) = " << ftXOptNM << "\n\n";

  printLine();
  std::cout << "\nsgpp::optimization example program terminated.\n";

  return 0;
}
