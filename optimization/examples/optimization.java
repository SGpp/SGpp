// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_optimization_java optimization.java
 *
 * On this page, we look at an example application of the sgpp::optimization module.
 * Versions of the example are given in all languages
 * currently supported by SG++: C++, Python, Java, and MATLAB.
 *
 * The example interpolates a bivariate test function like the \ref example_tutorial_cpp example.
 * However, we use B-splines here instead to obtain a smoother interpolant.
 * The resulting sparse grid function is then minimized with the method of steepest descent.
 * For comparison, we also minimize the objective function with Nelder-Mead's method.
 *
 * The example uses the following external class (stored in <tt>ExampleFunction.java</tt>).
 * The function \f$f\colon [0, 1]^d \to \mathbb{R}\f$ to be minimized
 * is called <i>objective function</i> and has to derive from sgpp.OptScalarFunction.
 * In the constructor, we give the dimensionality of the domain
 * (in this case \f$d = 2\f$).
 * The eval method evaluates the objective function and returns the function
 * value \f$f(\vec{x})\f$ for a given point \f$\vec{x} \in [0, 1]^d\f$.
 * \include ExampleFunction.java
 *
 * The actual example looks follows.
 * \dontinclude optimization.java
 */
public class optimization {
  static void printLine() {
    System.out.println("----------------------------------------" +
                       "----------------------------------------");
  }

  static String numToStr(double number) {
    return new java.text.DecimalFormat("#.#####").format(number);
  }

  public static void main(String[] args) {
    /**
     * At the beginning of the program, we have to load the shared library object file.
     * We can do so by using <tt>java.lang.System.load</tt> or
     * <tt>sgpp.LoadJSGPPLib.loadJSGPPLib</tt>.
     * Also, we disable OpenMP within jsgpp since it interferes with SWIG's director feature.
     */
    sgpp.LoadJSGPPLib.loadJSGPPLib();
    sgpp.jsgpp.omp_set_num_threads(1);

    System.out.println("sgpp::optimization example program started.\n");
    // increase verbosity of the output
    sgpp.OptPrinter.getInstance().setVerbosity(2);

    /**
     * Here, we define some parameters: objective function, dimensionality,
     * B-spline degree, maximal number of grid points, and adaptivity.
     */
    // objective function
    ExampleFunction f = new ExampleFunction();
    // dimension of domain
    final long d = f.getNumberOfParameters();
    // B-spline degree
    final long p = 3;
    // maximal number of grid points
    final long N = 30;
    // adaptivity of grid generation
    final double gamma = 0.95;

    /**
     * First, we define a grid with modified B-spline basis functions and
     * an iterative grid generator, which can generate the grid adaptively.
     */
    sgpp.Grid grid = sgpp.Grid.createModBsplineGrid(d, p);
    sgpp.OptIterativeGridGeneratorRitterNovak gridGen =
        new sgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma);

    /**
     * With the iterative grid generator, we generate adaptively a sparse grid.
     */
    printLine();
    System.out.println("Generating grid...\n");

    if (!gridGen.generate()) {
      System.out.println("Grid generation failed, exiting.");
      return;
    }

    /**
     * Then, we hierarchize the function values to get hierarchical B-spline
     * coefficients of the B-spline sparse grid interpolant
     * \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
     */
    printLine();
    System.out.println("Hierarchizing...\n");
    final sgpp.DataVector functionValues = gridGen.getFunctionValues();
    sgpp.DataVector coeffs = new sgpp.DataVector(functionValues.getSize());
    sgpp.OptHierarchisationSLE hierSLE = new sgpp.OptHierarchisationSLE(grid);
    sgpp.OptAutoSLESolver sleSolver = new sgpp.OptAutoSLESolver();

    // solve linear system
    if (!sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs)) {
      System.out.println("Solving failed, exiting.");
      return;
    }

    /**
     * We define the interpolant \f$\tilde{f}\f$ and its gradient
     * \f$\nabla\tilde{f}\f$ for use with the gradient method (steepest descent).
     * Of course, one can also use other optimization algorithms from
     * sgpp::optimization::optimizer.
     */
    printLine();
    System.out.println("Optimizing smooth interpolant...\n");
    sgpp.OptInterpolantScalarFunction ft =
        new sgpp.OptInterpolantScalarFunction(grid, coeffs);
    sgpp.OptInterpolantScalarFunctionGradient ftGradient =
        new sgpp.OptInterpolantScalarFunctionGradient(grid, coeffs);
    sgpp.OptGradientDescent gradientDescent =
        new sgpp.OptGradientDescent(ft, ftGradient);
    sgpp.DataVector x0 = new sgpp.DataVector(d);
    double fX0;
    double ftX0;

    /**
     * The gradient method needs a starting point.
     * We use a point of our adaptively generated sparse grid as starting point.
     * More specifically, we use the point with the smallest
     * (most promising) function value and save it in x0.
     */
    {
      sgpp.GridStorage gridStorage = grid.getStorage();

      // index of grid point with minimal function value
      int x0Index = 0;
      fX0 = functionValues.get(0);

      for (int i = 1; i < functionValues.getSize(); i++) {
        if (functionValues.get(i) < fX0) {
          fX0 = functionValues.get(i);
          x0Index = i;
        }
      }

      x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));
      ftX0 = ft.eval(x0);
    }

    System.out.println("x0 = " + x0);
    System.out.println("f(x0) = " + numToStr(fX0) +
                       ", ft(x0) = " + numToStr(ftX0) + "\n");

    /**
     * We apply the gradient method and print the results.
     */
    gradientDescent.setStartingPoint(x0);
    gradientDescent.optimize();
    sgpp.DataVector xOpt = gradientDescent.getOptimalPoint();
    final double ftXOpt = gradientDescent.getOptimalValue();
    final double fXOpt = f.eval(xOpt);

    System.out.println("\nxOpt = " + xOpt);
    System.out.println("f(xOpt) = " + numToStr(fXOpt) +
                       ", ft(xOpt) = " + numToStr(ftXOpt) + "\n");

    /**
     * For comparison, we apply the classical gradient-free Nelder-Mead method
     * directly to the objective function \f$f\f$.
     */
    printLine();
    System.out.println("Optimizing objective function (for comparison)...\n");

    sgpp.OptNelderMead nelderMead = new sgpp.OptNelderMead(f, 1000);
    nelderMead.optimize();
    sgpp.DataVector xOptNM = nelderMead.getOptimalPoint();
    final double fXOptNM = nelderMead.getOptimalValue();
    final double ftXOptNM = ft.eval(xOptNM);

    System.out.println("\nxOptNM = " + xOptNM);
    System.out.println("f(xOptNM) = " + numToStr(fXOptNM) +
                       ", ft(xOptNM) = " + numToStr(ftXOptNM) + "\n");

    printLine();
    System.out.println("\nsgpp::optimization example program terminated.");
  }
}  // end of class

/**
 * The example program outputs the following results:
 * \verbinclude optimization.output.txt
 *
 * We see that both the gradient-based optimization of the smooth sparse grid
 * interpolant and the gradient-free optimization of the objective function
 * find reasonable approximations of the minimum, which lies at
 * \f$(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)\f$.
 */
