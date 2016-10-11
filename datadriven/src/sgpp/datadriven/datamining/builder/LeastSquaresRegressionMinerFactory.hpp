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

#include "MinerFactory.hpp"

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class LeastSquaresRegressionMinerFactory : public MinerFactory {
 public:
  LeastSquaresRegressionMinerFactory();
  virtual ~LeastSquaresRegressionMinerFactory();

  virtual SparseGridMiner* buildMiner(const std::string& path);

 private:
  virtual DataSource* createDataSource(const DataMiningConfigParser& parser);
  virtual ModelFittingBase* createFitter(const DataMiningConfigParser& parser);
  virtual Scorer* createScorer(const DataMiningConfigParser& parser);
};

} /* namespace datadriven */
} /* namespace sgpp */
