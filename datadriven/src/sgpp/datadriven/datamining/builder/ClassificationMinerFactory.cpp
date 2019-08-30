/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ClassificationMinerFactory.cpp
 *
 *  Created on: Jul 2, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ClassificationFitterFactory.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase* ClassificationMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationClassification config{};
  config.readParams(parser);
  return new ModelFittingClassification(config);
}

FitterFactory *ClassificationMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  return new ClassificationFitterFactory(parser);
}
} /* namespace datadriven */
} /* namespace sgpp */
