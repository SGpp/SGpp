// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_fuzzy_cpp Fuzzy Extension Principle (C++)
 *
 * We consider an example for the application of the fuzzy extension principle to
 * fuzzy input uncertainties to obtain a fuzzy output interval using three different methods:
 * optimization of the objective function, optimization of a piecewise linear sparse grid
 * surrogate, and optimization of a surrogate with B-splines on sparse grids.
 * This example is available in C++ and in Python.
 */

/**
 * First, we include all the necessary headers, including those of the sgpp::base and
 * sgpp::optimization module.
 */
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>
#include <vector>

/**
 * Now, we can start with the \c main function.
 */
int main() {
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  /**
   * Here, we define some parameters: objective function, dimensionality,
   * B-spline degree, maximal number of grid points, and adaptivity.
   */
  // objective function
  sgpp::optimization::test_problems::Branin01Objective f;
  // dimension of domain
  const size_t d = f.getNumberOfParameters();
  // B-spline degree
  const size_t p = 3;
  // boundary parameter for the sparse grid
  const size_t b = 1;
  // level of the regular sparse grid
  const size_t n = 5;
  // accuracy of the extension principle
  const size_t numberOfAlphaSegments = 100;

  /**
   * We use regular sparse grids for the sparse grid surrogates.
   */
  std::unique_ptr<sgpp::base::Grid> gridBSpline;
  std::unique_ptr<sgpp::base::Grid> gridLinear;
  std::unique_ptr<sgpp::base::InterpolantScalarFunction> fInterpLinear;

  std::cout << "Constructing the sparse grids...\n";

  gridBSpline.reset(new sgpp::base::BsplineBoundaryGrid(d, p, b));
  gridBSpline->getGenerator().regular(n);

  gridLinear.reset(new sgpp::base::LinearBoundaryGrid(d, b));
  gridLinear->getGenerator().regular(n);

  const size_t N = gridBSpline->getSize();
  sgpp::base::GridStorage& gridStorage = gridBSpline->getStorage();

  sgpp::base::DataVector functionValues(N);
  sgpp::base::DataVector x(d);

  for (size_t k = 0; k < N; k++) {
    gridStorage[k].getStandardCoordinates(x);
    functionValues[k] = f.eval(x);
  }

  /**
   * For the hierarchization for the B-spline surrogate, we solve the corresponding
   * system of linear equations and create the interpolant and its gradient.
   */
  std::cout << "Hierarchizing (B-spline coefficients)...\n";

  std::unique_ptr<sgpp::base::InterpolantScalarFunction> fInterpBSpline;
  std::unique_ptr<sgpp::base::InterpolantScalarFunctionGradient> fInterpBSplineGradient;
  std::unique_ptr<sgpp::base::InterpolantScalarFunctionHessian> fInterpBSplineHessian;

  {
    sgpp::base::DataVector surpluses(N);
    sgpp::base::HierarchisationSLE hierSLE(*gridBSpline);
    sgpp::base::sle_solver::Auto sleSolver;

    if (!sleSolver.solve(hierSLE, functionValues, surpluses)) {
      std::cout << "Solving failed, exiting.\n";
      return 1;
    }

    fInterpBSpline.reset(
        new sgpp::base::InterpolantScalarFunction(*gridBSpline, surpluses));
    fInterpBSplineGradient.reset(
        new sgpp::base::InterpolantScalarFunctionGradient(*gridBSpline, surpluses));
  }

  /**
   * The piecewise linear interpolant is only continuous, but not continuously differentiable.
   * Therefore, we do not use the discontinuous gradient for gradient-based optimization,
   * but only use gradient-free optimization methods.
   */
  std::cout << "Hierarchizing (linear coefficients)...\n";

  {
    sgpp::base::DataVector surpluses(N);
    sgpp::base::HierarchisationSLE hierSLE(*gridLinear);
    sgpp::base::sle_solver::Auto sleSolver;

    if (!sleSolver.solve(hierSLE, functionValues, surpluses)) {
      std::cout << "Solving failed, exiting.\n";
      return 1;
    }

    fInterpLinear.reset(new sgpp::base::InterpolantScalarFunction(*gridLinear, surpluses));
  }

  std::cout << "\n";

  /**
   * Now we define the fuzzy input intervals.
   */
  sgpp::optimization::TriangularFuzzyInterval x0Fuzzy(0.25, 0.375, 0.125, 0.25);
  sgpp::optimization::QuasiGaussianFuzzyNumber x1Fuzzy(0.5, 0.125, 3.0);
  std::vector<sgpp::optimization::FuzzyInterval*> xFuzzy{&x0Fuzzy, &x1Fuzzy};

  /**
   * Finally, we can apply the fuzzy extension principle. First, we apply it directly
   * to the objective function to obtain a reference solution (fuzzy interval). Note that
   * usually the objective function is too expensive to use it directly in real-world
   * scenarios.
   */
  sgpp::optimization::optimizer::MultiStart optimizerExact(f, 10000, 100);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleExact(
      optimizerExact, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyExact(
      extensionPrincipleExact.apply(xFuzzy));
  std::cout << "L2 norm of exact solution: " << yFuzzyExact->computeL2Norm() << "\n";

  /**
   * For the piecewise linear and for the B-spline solution, we compute the relative
   * \f$L^2\f$ error to the exact solution computed previously.
   */
  sgpp::optimization::optimizer::MultiStart optimizerLinear(*fInterpLinear, 10000, 100);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleLinear(
      optimizerLinear, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyLinear(
      extensionPrincipleLinear.apply(xFuzzy));
  std::cout << "Relative L2 error of piecewise linear solution: " <<
      yFuzzyExact->computeRelativeL2Error(*yFuzzyLinear) << "\n";

  /**
   * For B-splines, we use gradient descent as our optimization method.
   */
  sgpp::optimization::optimizer::AdaptiveGradientDescent localOptimizer(
      *fInterpBSpline, *fInterpBSplineGradient);
  sgpp::optimization::optimizer::MultiStart optimizerBSpline(localOptimizer);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleBSpline(
      optimizerBSpline, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyBSpline(
      extensionPrincipleBSpline.apply(xFuzzy));
  std::cout << "Relative L2 error of B-spline solution: " <<
      yFuzzyExact->computeRelativeL2Error(*yFuzzyBSpline) << "\n";

  return 0;
}

/**
 * The example outputs something similar to the following:
 * \verbinclude fuzzy.output.txt
 *
 * The exact output is not deterministic (despite setting the seed of the
 * random number generator at the start of <tt>main</tt>), since the \c MultiStart
 * optimizer calls are executed in parallel for the different confidence intervals.
 *
 * We see that the relative \f$L^2\f$ error is over three orders of magnitude smaller
 * for the B-spline solution compared to the piecewise linear solution.
 * This is due to the higher approximation quality of the B-splines and
 * to the gradient-based optimization.
 */
