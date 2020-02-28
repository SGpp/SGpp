// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Concrete Factory that builds an instance of #sgpp::datadriven::SparseGridMiner for classification
 */
class ClassificationMinerFactory : public MinerFactory {
 public:
  /**
   * Default constructor
   */
  ClassificationMinerFactory() = default;

 private:
  /**
   * Build an instance of a #sgpp::datadriven::ModelFittingBase object as specified in the
   * configuration
   * file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured fitter (instance of a #sgpp::datadriven::ModelFittingBase object) as
   * specified in the
   * configuration file.
   */
  ModelFittingBase* createFitter(const DataMiningConfigParser& parser) const override;

  FitterFactory* createFitterFactory(const DataMiningConfigParser& parser) const override;

  /* Factory method to build a visualizer instance base on a configuration file.
   * @param parser the datamining configuration parser instance to create the scorer from
   * @return the scorer instance
   */
  Visualizer* createVisualizer(const DataMiningConfigParser& parser) const override;
};

} /* namespace datadriven */
} /* namespace sgpp */
