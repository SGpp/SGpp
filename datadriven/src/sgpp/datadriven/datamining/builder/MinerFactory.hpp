/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * MinerFactory.hpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>

#include <string>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Abstract factory to build different kinds of Miners based on a configuration which is parsed from
 * a file. A miner consists of a data source, a fitter and a scorer. A concrete Factory class has to
 * implement the required interfaces.
 */
class MinerFactory {
 public:
  /**
   * Default constructor
   */
  MinerFactory() = default;

  /**
   * Virtual destructor
   */
  virtual ~MinerFactory() = default;

  /**
   * Factory method to build a miner object based on a configuration file.
   * @param path Path to a configuration file that defines the structure of the miner object.
   */
  virtual SparseGridMiner* buildMiner(const std::string& path) const = 0;

  virtual HyperparameterOptimizer* buildHPO(const std::string& path) const = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
