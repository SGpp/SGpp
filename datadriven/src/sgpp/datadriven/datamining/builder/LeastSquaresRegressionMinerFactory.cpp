/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFactory.cpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase* LeastSquaresRegressionMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  return new ModelFittingLeastSquares(config);
}

FitterFactory *LeastSquaresRegressionMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  return new LeastSquaresRegressionFitterFactory(parser);
}
} /* namespace datadriven */
} /* namespace sgpp */
