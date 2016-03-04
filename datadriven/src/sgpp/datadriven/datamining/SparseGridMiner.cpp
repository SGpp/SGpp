/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: franzefn
 */

#include <sgpp/datadriven/datamining/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/ARFFWrapper.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/SimpleSplittingScorer.hpp>
#include <sgpp/datadriven/datamining/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/MSE.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

SparseGridMiner::SparseGridMiner(datadriven::DataMiningConfiguration pconfig)
    : scorer(NULL), config(pconfig) {}

SparseGridMiner::~SparseGridMiner() {}

void SparseGridMiner::run() {
  // 1. find lambda
  double threshold = config["scorer_threshold"].getDouble();
  uint64_t maxRefinenum = config["maxRefinenum"].getUInt();

  // 2. read data set
  datadriven::ARFFWrapper dataWrapper(config);
  datadriven::Dataset dataset = dataWrapper.allSamples();

  // 3. generate learner

  std::shared_ptr<datadriven::ModelFittingDensityEstimation> fitter =
      std::make_shared<datadriven::ModelFittingDensityEstimation>(config);

  // 4. load simple splitting scorer
  std::shared_ptr<datadriven::Metric> metric = std::make_shared<datadriven::MSE>();
  datadriven::SimpleSplittingScorer scorer(metric, fitter, config);

  double score = scorer.getScore(dataset);
  uint32_t iter = 0;
  while (score > threshold && iter < maxRefinenum) {
    fitter->refine();
    score = scorer.getScore(dataset);
    iter++;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
