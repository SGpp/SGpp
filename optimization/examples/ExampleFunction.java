// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

public class ExampleFunction extends sgpp.OptObjectiveFunction {
  public ExampleFunction() {
    super(2);
  }

  public double eval(sgpp.DataVector x) {
    if ((x.get(0) >= 0.0) && (x.get(0) <= 1.0) &&
        (x.get(1) >= 0.0) && (x.get(1) <= 1.0)) {
      return Math.sin(8.0 * x.get(0)) + Math.sin(7.0 * x.get(1));
    } else {
      return Double.POSITIVE_INFINITY;
    }
  }
}
