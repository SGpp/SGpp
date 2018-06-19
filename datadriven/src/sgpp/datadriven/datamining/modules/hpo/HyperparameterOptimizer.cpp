/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HyperparameterOptimizer.cpp
 *
 * Created on: Jan 22, 2018
 *     Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>

#include <sgpp/datadriven/datamining/modules/hpo/HPOScorerFactory.hpp>

#include <vector>
#include <string>
#include <limits>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource *dataSource,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
    : fitterFactory(fitterFactory) {
  HPOScorerFactory scorerFactory;
  hpoScorer.reset(dynamic_cast<HPOScorer *>(scorerFactory.buildScorer(parser)));
  config.setupDefaults();
  parser.getHPOConfig(config);
  std::unique_ptr<DataSource> ds(dataSource);
  trainData.reset(ds->getNextSamples());
  if (!parser.hasScorerTestset()) {
    trainData.reset(hpoScorer->prepareTestData(*trainData));
  }
  if (config.getNTrainSamples() > 0
      && static_cast<size_t>(config.getNTrainSamples()) <= trainData->getNumberInstances()) {
    Dataset *resize =
        new Dataset(static_cast<size_t>(config.getNTrainSamples()), trainData->getDimension());
    hpoScorer->resizeTrainData(*trainData, *resize);
    trainData.reset(resize);
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
