/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * MinerPrototype.cpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::LeastSquaresRegressionMinerFactory;
using sgpp::datadriven::SparseGridMiner;

int main(int argc, char **argv) {
  // read input
  std::string path;
  if (argc != 2) {
    std::cout << "No or bad path given, aborting" << std::endl;
    exit(-1);
  } else {
    path = std::string(argv[1]);
  }

  LeastSquaresRegressionMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(path));
  miner->learn();
}
