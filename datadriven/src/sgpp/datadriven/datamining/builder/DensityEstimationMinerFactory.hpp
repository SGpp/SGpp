/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DensityEstimationMinerFactory.hpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#pragma once

#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Concrete Factory that builds an instance of #sgpp::datadriven::SparseGridMiner for Density Estimation
 */
class DensityEstimationMinerFactory : public MinerFactory {
 public:
  /**
   * Default constructor
   */
  DensityEstimationMinerFactory() = default;

  /**
   * Build an instance of #sgpp::datadriven::SparseGridMiner for Density Estimation based on
   * specification from a configuration file.
   * @param path Path to a configuration file that defines the structure of the miner object.
   */
  virtual SparseGridMiner* buildMiner(const std::string& path) const;

 private:
  /**
   * Build an instance of a #sgpp::datadriven::DataSource object as specified in the configuration
   * file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured instance of a #sgpp::datadriven::DataSource object as specified in the
   * configuration file.
   */
  virtual DataSource* createDataSource(const DataMiningConfigParser& parser) const;

  /**
   * Build an instance of a #sgpp::datadriven::ModelFittingBase object as specified in the
   * configuration
   * file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured fitter (instance of a #sgpp::datadriven::ModelFittingBase object) as
   * specified in the
   * configuration file.
   */
  virtual ModelFittingBase* createFitter(const DataMiningConfigParser& parser) const;

  /**
   * Build an instance of a #sgpp::datadriven::Scorer object as specified in the configuration
   * file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured instance of a #sgpp::datadriven::Scorer object as specified in the
   * configuration file.
   */
  virtual Scorer* createScorer(const DataMiningConfigParser& parser) const;
};

} /* namespace datadriven */
} /* namespace sgpp */
