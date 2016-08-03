// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
* \page example_dataVectorSerializeDemo_cpp dataVectorSerializeDemo.cpp
*
* This example shows how to initialize, serialize, and deserialize a DataVector object.
*/

#include <iostream>

#include "sgpp/base/datatypes/DataVector.hpp"

using sgpp::base::DataVector;

int main() {
  /**
   * We create a DataVector with 3 elements and initialize the values.
   */
  DataVector v;
  v.append(1.0);
  v.append(2.0);
  v.append(3.0);

  /**
   * Store it to a file using the toFile() method
   */
  v.toFile("dataVectorTest.vec");

  /**
   * Restore it from a file using the fromFile() class method.
   */
  DataVector w = DataVector::fromFile("dataVectorTest.vec");

  /**
   * Store again (for no reason).
   */
  w.toFile("dataVectorTest2.vec");
}
