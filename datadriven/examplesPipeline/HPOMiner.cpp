// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>


#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::UniversalMinerFactory;
using sgpp::datadriven::SparseGridMiner;

/**
 * This example demonstrates how a #sgpp::datadriven::SparseGridMiner is constructed using a
 * configuration file and how it is then used.
 */
int main(int argc, char **argv) {
  /**
   * use immediately invoked lambda expression to get the path to a configuration file.
   */
  const std::string path = [argc, &argv]() {
    if (argc != 2) {
      std::cout << "No or bad path given, aborting\n";
      exit(1);
    } else {
      return std::string{argv[1]};
    }
  }();

  /**
   * We need a factory class to actually build the #sgpp::datadriven::SparseGridMiner.
   */
  UniversalMinerFactory factory;
  /**
   * The miner object is constructed by the factory from a supplied configuration file.
   */
  auto hpo = std::unique_ptr<sgpp::datadriven::HyperparameterOptimizer>(factory.buildHPO(path));
  /**
   * Once we have a configured miner object, we can start the learning process.
   */
  hpo->run(true);
}
