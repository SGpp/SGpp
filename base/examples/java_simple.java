// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// import all packages
//import sgpp.*;

// or, better, include only the ones needed
import sgpp.LoadJSGPPLib;
import sgpp.DataVector;
import sgpp.GridGenerator;
import sgpp.GridStorage;
import sgpp.Grid;
import sgpp.GridIndex;
import sgpp.jsgpp;
import sgpp.OperationEval;

public class java_simple {
  private static double f(double x0, double x1) {
    return 16.0*(x0-1.0)*x0 * (x1-1.0)*x1;
  }

  public static void main(String[] args) {
    // Two working possiblities for loading the shared object file:
    //java.lang.System.load("/path/to/SGpp/trunk/lib/jsgpp/libjsgpp.so");
    sgpp.LoadJSGPPLib.loadJSGPPLib();

    // create a two-dimensional piecewise bilinear grid
    int dim = 2;
    Grid grid = Grid.createLinearGrid(dim);
    GridStorage gridStorage = grid.getStorage();
    System.out.println("dimensionality:         " + gridStorage.dim());

    // create regular grid, level 3
    int level = 3;
    GridGenerator gridGen = grid.createGridGenerator();
    gridGen.regular(level);

    System.out.println("number of grid points:  " + gridStorage.size());

    // create coefficient vector
    DataVector alpha = new DataVector(gridStorage.size());
    alpha.setAll(0.0);
    System.out.println("length of alpha vector: " + alpha.getSize());

    // set function values in alpha
    GridIndex gp = new GridIndex();

    for (int i = 0; i < gridStorage.size(); i++) {
      gp = gridStorage.get(i);
      alpha.set(i, f(gp.getCoord(0), gp.getCoord(1)));
    }
    System.out.println("alpha before hierarchization: " + alpha.toString());

    // hierarchize
    jsgpp.createOperationHierarchisation(grid).doHierarchisation(alpha);
    System.out.println("alpha after hierarchization:  " + alpha.toString());

    // evaluate
    DataVector p = new DataVector(dim);
    p.set(0, 0.52);
    p.set(1, 0.73);
    OperationEval opEval = jsgpp.createOperationEval(grid);

    System.out.println("u(0.52, 0.73) = " + opEval.eval(alpha,p));
  }
}
