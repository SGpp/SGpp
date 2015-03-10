// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

public class java_example {
  static void printLine() {
    System.out.println("----------------------------------------" +
                       "----------------------------------------");
  }

  static String numToStr(double number) {
    return new java.text.DecimalFormat("#.#####").format(number);
  }

  public static void main(String[] args) {
    // load jsgpp
    sgpp.LoadJSGPPLib.loadJSGPPLib();
    // disable OpenMP multi-threading within jsgpp
    // (interferes with SWIG's director feature)
    sgpp.jsgpp.omp_set_num_threads(1);
    // increase output verbosity
    sgpp.jsgpp.getOptPrinterInstance().setVerbosity(2);

    System.out.println("SGPP::optimization example program started.\n");

    // objective function
    ExampleFunction f = new ExampleFunction();
    // dimension of domain
    final long d = f.getDimension();
    // B-spline degree
    final long p = 3;
    // maximal number of grid points
    final long N = 30;
    // adaptivity of grid generation
    final double gamma = 0.95;

    sgpp.Grid grid = sgpp.Grid.createModBsplineGrid(f.getDimension(), p);
    sgpp.OptIterativeGridGeneratorRitterNovak gridGen =
        new sgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma);

    // ////////////////////////////////////////////////////////////////////////
    // GRID GENERATION
    // ////////////////////////////////////////////////////////////////////////

    printLine();
    System.out.println("Generating grid...\n");

    if (!gridGen.generate()) {
      System.out.println("Grid generation failed, exiting.");
      return;
    }

    // ////////////////////////////////////////////////////////////////////////
    // HIERARCHIZATION
    // ////////////////////////////////////////////////////////////////////////

    printLine();
    System.out.println("Hierarchizing...\n");
    sgpp.DoubleVector coeffs = new sgpp.DoubleVector();
    sgpp.OptHierarchisationSLE hierSLE = new sgpp.OptHierarchisationSLE(grid);
    sgpp.OptAutoSLESolver sleSolver = new sgpp.OptAutoSLESolver();

    // solve linear system
    if (!sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs)) {
      System.out.println("Solving failed, exiting.");
      return;
    }

    // convert std::vector to sgpp.DataVector
    sgpp.DataVector coeffsDV = new sgpp.DataVector(coeffs);

    // ////////////////////////////////////////////////////////////////////////
    // OPTIMIZATION OF THE SMOOTH INTERPOLANT
    // ////////////////////////////////////////////////////////////////////////

    printLine();
    System.out.println("Optimizing smooth interpolant...\n");
    sgpp.OptInterpolantFunction ft =
        new sgpp.OptInterpolantFunction(d, grid, coeffsDV);
    sgpp.OptInterpolantGradient ftGradient =
        new sgpp.OptInterpolantGradient(d, grid, coeffsDV);
    sgpp.OptGradientMethod gradientMethod =
        new sgpp.OptGradientMethod(ft, ftGradient);
    sgpp.DoubleVector x0 = new sgpp.DoubleVector(d);
    double fX0;
    double ftX0;

    // determine best grid point as starting point
    {
      final sgpp.DoubleVector functionValues = gridGen.getFunctionValues();
      sgpp.GridStorage gridStorage = grid.getStorage();

      // index of grid point with minimal function value
      int x0Index = 0;
      fX0 = functionValues.get(0);

      for (int i = 1; i < functionValues.size(); i++) {
        if (functionValues.get(i) < fX0) {
          fX0 = functionValues.get(i);
          x0Index = i;
        }
      }

      for (int t = 0; t < d; t++) {
        x0.set(t, gridStorage.get(x0Index).getCoord(t));
      }

      ftX0 = ft.eval(x0);
    }

    System.out.println("x0 = " + new sgpp.DataVector(x0));
    System.out.println("f(x0) = " + numToStr(fX0) +
                       ", ft(x0) = " + numToStr(ftX0) + "\n");

    gradientMethod.setStartingPoint(x0);
    sgpp.DoubleVector xOpt = new sgpp.DoubleVector();
    final double ftXOpt = gradientMethod.optimize(xOpt);
    final double fXOpt = f.eval(xOpt);

    System.out.println("\nxOpt = " + new sgpp.DataVector(xOpt));
    System.out.println("f(xOpt) = " + numToStr(fXOpt) +
                       ", ft(xOpt) = " + numToStr(ftXOpt) + "\n");

    // ////////////////////////////////////////////////////////////////////////
    // NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
    // ////////////////////////////////////////////////////////////////////////

    printLine();
    System.out.println("Optimizing objective function (for comparison)...\n");

    sgpp.OptNelderMead nelderMead = new sgpp.OptNelderMead(f, 1000);
    sgpp.DoubleVector xOptNM = new sgpp.DoubleVector();
    final double fXOptNM = nelderMead.optimize(xOptNM);
    final double ftXOptNM = ft.eval(xOptNM);

    System.out.println("\nxOptNM = " + new sgpp.DataVector(xOptNM));
    System.out.println("f(xOptNM) = " + numToStr(fXOptNM) +
                       ", ft(xOptNM) = " + numToStr(ftXOptNM) + "\n");

    printLine();
  
    System.out.println("\nSGPP::optimization example program terminated.");
  }
}
