// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
   * \page example_dataVectorSerializeDemo_cpp Using the DataVector object
   *
   * This example shows how to initialize a DataVector object, store it to a file
   * and then to restore it back.
   */

#include <sgpp/base/datatypes/DataVector.hpp>

#include <iostream>


using sgpp::base::DataVector;

int main() {
  /**
   * We create a DataVector and fill it with 3 values
   */
  DataVector v;
  v.append(1.0);
  v.append(2.0);
  v.append(3.0);

  /**
   * We save the DataVector to file
   */
  v.toFile("dataVectorTest.vec");

  /**
   * We load the saved DataVector and then save it again (for no particular reason)
   */
  DataVector w = DataVector::fromFile("dataVectorTest.vec");
  w.toFile("dataVectorTest2.vec");
}
