// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_optimization_cpp optimization.cpp
 *
 * On this page, we look at an example application of the sgpp::optimization module.
 * Versions of the example are given in all languages
 * currently supported by SG++: C++, Python, Java, and MATLAB.
 *
 * The example interpolates a bivariate test function like the \ref example_tutorial_cpp example.
 * However, we use B-splines here instead to obtain a smoother interpolant.
 * The resulting sparse grid function is then minimized with the method of steepest descent.
 * For comparison, we also minimize the objective function with Nelder-Mead's method.
 */

/**
 * First, we include all the necessary headers, including those of the sgpp::base and
 * sgpp::optimization module.
 */
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>

/**
 * The function \f$f\colon [0, 1]^d \to \mathbb{R}\f$ to be minimized
 * is called <i>objective function</i> and has to derive from
 * sgpp::optimization::ScalarFunction.
 * In the constructor, we give the dimensionality of the domain
 * (in this case \f$d = 2\f$).
 * The eval method evaluates the objective function and returns the function
 * value \f$f(\vec{x})\f$ for a given point \f$\vec{x} \in [0, 1]^d\f$.
 * The clone method returns a std::unique_ptr to a clone of the object
 * and is used for parallelization (in case eval is not thread-safe).
 */
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

/**
 * Now, we can start with the \c main function.
 */
void printLine() {
  std::cout << "----------------------------------------"
               "----------------------------------------\n";
}

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "sgpp::optimization example program started.\n\n";
  // increase verbosity of the output
  sgpp::optimization::Printer::getInstance().setVerbosity(2);

  /**
   * Here, we define some parameters: objective function, dimensionality,
   * B-spline degree, maximal number of grid points, and adaptivity.
   */
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

  /**
   * First, we define a grid with modified B-spline basis functions and
   * an iterative grid generator, which can generate the grid adaptively.
   */
  sgpp::base::ModBsplineGrid grid(d, p);
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(f, grid, N, gamma);

  /**
   * With the iterative grid generator, we generate adaptively a sparse grid.
   */
  printLine();
  std::cout << "Generating grid...\n\n";

  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  /**
   * Then, we hierarchize the function values to get hierarchical B-spline
   * coefficients of the B-spline sparse grid interpolant
   * \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
   */
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

  /**
   * We define the interpolant \f$\tilde{f}\f$ and its gradient
   * \f$\nabla\tilde{f}\f$ for use with the gradient method (steepest descent).
   * Of course, one can also use other optimization algorithms from
   * sgpp::optimization::optimizer.
   */
  printLine();
  std::cout << "Optimizing smooth interpolant...\n\n";
  sgpp::optimization::InterpolantScalarFunction ft(grid, coeffs);
  sgpp::optimization::InterpolantScalarFunctionGradient ftGradient(grid, coeffs);
  sgpp::optimization::optimizer::GradientDescent gradientDescent(ft, ftGradient);
  sgpp::base::DataVector x0(d);
  double fX0;
  double ftX0;

  /**
   * The gradient method needs a starting point.
   * We use a point of our adaptively generated sparse grid as starting point.
   * More specifically, we use the point with the smallest
   * (most promising) function value and save it in x0.
   */
  {
    sgpp::base::GridStorage& gridStorage = grid.getStorage();

    // index of grid point with minimal function value
    size_t x0Index =
        std::distance(functionValues.getPointer(),
                      std::min_element(functionValues.getPointer(),
                                       functionValues.getPointer() + functionValues.getSize()));

    x0 = gridStorage.getCoordinates(gridStorage[x0Index]);

    fX0 = functionValues[x0Index];
    ftX0 = ft.eval(x0);
  }

  std::cout << "x0 = " << x0.toString() << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  /**
   * We apply the gradient method and print the results.
   */
  gradientDescent.setStartingPoint(x0);
  gradientDescent.optimize();
  const sgpp::base::DataVector& xOpt = gradientDescent.getOptimalPoint();
  const double ftXOpt = gradientDescent.getOptimalValue();
  const double fXOpt = f.eval(xOpt);

  std::cout << "\nxOpt = " << xOpt.toString() << "\n";
  std::cout << "f(xOpt) = " << fXOpt << ", ft(xOpt) = " << ftXOpt << "\n\n";

  /**
   * For comparison, we apply the classical gradient-free Nelder-Mead method
   * directly to the objective function \f$f\f$.
   */
  printLine();
  std::cout << "Optimizing objective function (for comparison)...\n\n";

  sgpp::optimization::optimizer::NelderMead nelderMead(f, 1000);
  nelderMead.optimize();
  sgpp::base::DataVector xOptNM = nelderMead.getOptimalPoint();
  const double fXOptNM = nelderMead.getOptimalValue();
  const double ftXOptNM = ft.eval(xOptNM);

  std::cout << "\nnxOptNM = " << xOptNM.toString() << "\n";
  std::cout << "f(xOptNM) = " << fXOptNM << ", ft(xOptNM) = " << ftXOptNM << "\n\n";

  printLine();
  std::cout << "\nsgpp::optimization example program terminated.\n";

  return 0;
}

/**
 * The example program outputs the following results:
 * \verbinclude optimization.output.txt
 *
 * We see that both the gradient-based optimization of the smooth sparse grid
 * interpolant and the gradient-free optimization of the objective function
 * find reasonable approximations of the minimum, which lies at
 * \f$(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)\f$.
 */
