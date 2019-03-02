/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SGDEMinerFromConfigFile.cpp
 *
 *  Created on: Jun 9, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;

/**
 * This example demonstrates how a #sgpp::datadriven::SparseGridMiner is constructed using a
 * configuration file and how it is then used.
 */
/*int main(int argc, char **argv) {*/
    /**
   * use immediately invoked lambda expression to get the path to a configuration file.
   */
/*  const std::string path = [argc, &argv]() {
#    if (argc != 2) {
#      std::cout << "No or bad path given, aborting\n";
#      exit(1);
#      return std::string{};
#    } else {
#      return std::string{argv[1]};
#    }
#  }();*/
int main(int argc) {
  const std::string path = "classificationMinerConfigCoarsening.json";

  /**
   * We need a factory class to actually build the #sgpp::datadriven::SparseGridMiner.
   */
  ClassificationMinerFactory factory;

  /**
   * The miner object is constructed by the factory from a supplied configuration file.
   */
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(path));
  /**
   * Once we have a configured miner object, we can start the learning process.
   */
  miner->learn(true);
  std::cout << std::endl;
}




