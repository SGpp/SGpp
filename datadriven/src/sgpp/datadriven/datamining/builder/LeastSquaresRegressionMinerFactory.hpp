/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFactory.hpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Concrete Factory that builds an instance of #sgpp::datadriven::SparseGridMiner for Least Squares
 * Regression
 */
class LeastSquaresRegressionMinerFactory : public MinerFactory {
 public:
  /**
   * Default constructor
   */
  LeastSquaresRegressionMinerFactory() = default;

  sgpp::datadriven::HyperparameterOptimizer *buildHPO(const std::string &path) const override;

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
};
} /* namespace datadriven */
} /* namespace sgpp */
