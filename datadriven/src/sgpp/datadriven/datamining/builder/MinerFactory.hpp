// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>

#include <string>
#include <vector>

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
  virtual SparseGridMiner* buildMiner(const std::string& path) const;

  virtual sgpp::datadriven::HyperparameterOptimizer* buildHPO(const std::string& path) const;

 protected:
  /**
   * Factory method to build a splitting based data source, i.e. a data source that splits data into
   * validation and training data.
   * @param parser the datamining configuration parser instance to create the data source from
   * @return the data source instance
   */
  virtual DataSourceSplitting* createDataSourceSplitting(
      const DataMiningConfigParser& parser) const;

  virtual std::vector<DataSourceSplitting*> createDataSourceSplittingTwoDatasets(
      const DataMiningConfigParser& parser) const {
    throw base::application_exception("This miner only allows a single dataSource instance");
  }

  /**
   * Factory method to build a cross validation data source, i.e. a data source that can separate
   * one fold from the data as validation set and use the rest for training
   * @param parser the datamining configuration parser instance to create the data source from
   * @return the data source instance
   */
  virtual DataSourceCrossValidation* createDataSourceCrossValidation(
      const DataMiningConfigParser& parser) const;

  /**
   * Build an instance of a #sgpp::datadriven::ModelFittingBase object as specified in the
   * configuration file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured fitter (instance of a #sgpp::datadriven::ModelFittingBase object) as
   * specified in the configuration file.
   */
  virtual ModelFittingBase* createFitter(const DataMiningConfigParser& parser) const = 0;

  virtual FitterFactory* createFitterFactory(const DataMiningConfigParser& parser) const = 0;

  /**
   * Factory method to build a scorer instance base d on a configuration file.
   * @param parser the datamining configuration parser instance to create the scorer from
   * @return the scorer instance
   */
  virtual Scorer* createScorer(const DataMiningConfigParser& parser) const;

  /**
   * Factory method to build a visualizer instance base on a configuration file.
   * @param parser the datamining configuration parser instance to create the scorer from
   * @return the scorer instance
   */
  virtual Visualizer* createVisualizer(const DataMiningConfigParser& parser) const = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
