/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ParallelClassificationMinerFromConfigFile.cpp
 *
 * Created on: Apr 02, 2019
 *     Author: Jan Schopohl
 */

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;

/**
 * This example demonstrates how a parallel (distributed) sgpp::datadriven::SparseGridMiner is
 * constructed using a configuration file and how it is then used. Based on
 * ClassificationMinerFromConfigFile.cpp
 */
int main(int argc, char **argv) {
  /**
   * Initialize MPI, BLACS and ScaLAPACK
   */
  sgpp::datadriven::BlacsProcessGrid::initializeBlacs();
  {
    /**
     * use immediately invoked lambda expression to get the path to a configuration file.
     */
    const std::string path = [argc, &argv]() {
      if (argc != 2) {
        std::cout << "No or bad path given, aborting\n";
        exit(1);
        return std::string{};
      } else {
        return std::string{argv[1]};
      }
    }();

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
  /**
   * exit BLACS, MPI (in outer scope to ensure all blacs grids were destructed before)
   */
  sgpp::datadriven::BlacsProcessGrid::exitBlacs();
}
