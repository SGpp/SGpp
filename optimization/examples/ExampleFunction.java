// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

public class ExampleFunction extends sgpp.OptObjectiveFunction {
  public ExampleFunction() {
    super(2);
  }

  public double eval(sgpp.DoubleVector x) {
    return Math.sin(8.0 * x.get(0)) + Math.sin(7.0 * x.get(1));
  }
}
