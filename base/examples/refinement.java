// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
  * \page example_refinement_java Refinement Example
  *
  * Here we demonstrate how to refine a grid. As a refinement indicator, we take the surpluses
  * of the grid points directly.
  * We start with a regular sparse grid of level 3 with linear basis functions and refine five
  * times. In each refinement step, we refine the grid point with the highest absolute surplus.
  *
  * The following example interpolates the (non-symmetric) function
  * f : [0,1]x[0,1] -> R, f(x,y) := 16 (x-1) * x * (y-1) * y
  *
  * The number of grid points is printed in each iteration.
  * After refinement, the surplusses have to be set for all new grid
  * points, i.e., the alpha-Vector has to be extended.
  *
  * For instructions on how to run the example, please see \ref installation.
  */

import sgpp.*;

class refinement {
  /**
    * function to interpolate. This is a two-dimensional parabola. - nonsymmetric(!)
    */
  private static double f(double x0, double x1) {
    return 16.0 * (x0 - 1.0) * x0 * (x1 - 1.0) * x1;
  }

  public static void main(String[] args) {
    sgpp.LoadJSGPPLib.loadJSGPPLib();

    /**
      * create a two-dimensional piecewise bi-linear grid
      */
    int dim = 2;
    Grid grid = Grid.createLinearGrid(dim);
    GridStorage gridStorage = grid.getStorage();
    System.out.println("dimensionality:                   " + gridStorage.getDimension());


    /**
      * create regular sparse grid, level 3
      */
    int level = 3;
    grid.getGenerator().regular(level);
    System.out.println("number of initial grid points:    " + gridStorage.getSize());


    /**
      * create coefficient vector with size corresponding to the grid size.
      * Initially, all the values are set to zero.
      */
    DataVector alpha = new DataVector(gridStorage.getSize());
    alpha.setAll(0.0);
    System.out.println("length of alpha vector:           " + alpha.getSize());

    /**
      * Obtain function values and refine adaptively 5 times
      */
    for (int step = 0; step < 5; step++) {
      // set function values in alpha
      for (int i = 0; i < gridStorage.getSize(); i++) {
        GridPoint gp = gridStorage.getPoint(i);
        alpha.set(i, f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)));
      }

      /**
        * Each time, we have to hierarchize the grid again, because in the previous interation,
        * new grid points have been added.
        */
      jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);

      /**
        * Refine a single grid point each time.
        * The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
        * Refining the point means, that all children of this point (if not already present) are
        * added to the grid. Also all missing parents are added (recursively).
        */
      SurplusRefinementFunctor functor = new SurplusRefinementFunctor(alpha, 1);
      grid.getGenerator().refine(functor);
      System.out.println("refinement step " + (step + 1) + ", new grid size: " + alpha.getSize());

      /**
        * Extend alpha vector (new entries uninitialized). Note that right now, the surplus vector
        * has the correct size again, but the values of the new points are set to zero. The correct
        * surplus values will be inserted in the next iteration of the refinement loop.
        */
      alpha.resize(gridStorage.getSize());
    }
  }
}

/**
 * This results in the following output:
 * \verbinclude refinement.output.txt
 *
 * There are clearly more efficient approaches than to set the function
 * values for all grid points and to hierarchize the whole grid each
 * time. But this works even where no efficient alternatives are
 * available and suffices for demonstration purposes.
 *
 * This use of the SurplusRefinementFunctor takes as arguments the
 * coefficient vector (it doesn't have to be the coefficient vector, it
 * could be something modified!) and the number of grid points to refine
 * (if available). It bases its refinement decision on the absolute
 * values of the vector's entries, choosing the largest ones. Other
 * refinement functors are available or can be implemented.
 */
