// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
   * \page example_dataMatrixSerializeDemo_cpp Using the DataMatrix object
   *
   * This example shows how to initialize a DataMatrix object, store it to a file
   * and then to restore it back.
   */

#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <iostream>


using sgpp::base::DataMatrix;

int main() {
  /**
     * We create a 2-by-2 matrix and fill it with values using the set() function.
     */
  DataMatrix m(2, 2);
  m.set(0, 0, 1.0);
  m.set(0, 1, 2.0);
  m.set(1, 0, 3.0);
  m.set(1, 1, 4.0);

  /**
     * Now we store the matrix to a file
     */
  m.toFile("dataMatrixTest.mat");

  /**
     * We load a DataMatrix from a file and store it again (for no particular reason).
     */
  DataMatrix m2 = DataMatrix::fromFile("dataMatrixTest.mat");

  m2.toFile("dataMatrixTest2.mat");
}
