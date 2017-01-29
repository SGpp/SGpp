/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * IChol.cpp
 *
 *  Created on: Jan 29, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/IChol.hpp>

#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::IChol;

int main(int argc, char** argv) {
  const auto sweeps = [argc, argv]() {
    if (argc != 2) {
      printf("invalid amount of sweeps specified");
      exit(1);
    } else {
      try {
        return std::stoul(argv[1]);
      } catch (std::invalid_argument& e) {
        printf("amount of sweeps could not be parsed.\n%s\n", e.what());
        exit(1);
      }
    }
  }();

  auto size = 2u;

  // initialize
  double a_val[]{8.0, 4.0, 4.0, 9.0};
  double b_val[]{2.0 * sqrt(2), 4.0, sqrt(2), sqrt(7)};

  DataVector norm{size};
  DataMatrix A{a_val, size, size};
  DataMatrix B{b_val, size, size};

  // decomp:

  IChol::normToUnitDiagonal(A, norm);
  IChol::decompose(A, sweeps);
  IChol::reaplyDiagonal(A, norm);

  printf("Chol\n\n%s", B.toString().c_str());

  printf("IChol\n\n%s\n\n", A.toString().c_str());
}
