// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>

#include "sgpp/base/datatypes/DataVector.hpp"

using sgpp::base::DataVector;

int main() {
  DataVector v;
  v.append(1.0);
  v.append(2.0);
  v.append(3.0);

  v.toFile("dataVectorTest.vec");

  DataVector w = DataVector::fromFile("dataVectorTest.vec");

  w.toFile("dataVectorTest2.vec");
}
