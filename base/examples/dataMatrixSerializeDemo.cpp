// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>

#include "sgpp/base/datatypes/DataMatrix.hpp"

using sgpp::base::DataMatrix;

int main() {
  DataMatrix m(2, 2);
  m.set(0, 0, 1.0);
  m.set(0, 1, 2.0);
  m.set(1, 0, 3.0);
  m.set(1, 1, 4.0);

  m.toFile("dataMatrixTest.mat");

  DataMatrix m2 = DataMatrix::fromFile("dataMatrixTest.mat");

  m2.toFile("dataMatrixTest2.mat");
}
