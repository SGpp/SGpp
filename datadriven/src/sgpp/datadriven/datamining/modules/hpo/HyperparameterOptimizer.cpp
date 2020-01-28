// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>


#include <vector>
#include <string>
#include <limits>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(SparseGridMiner* miner,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
    : miner(miner), fitterFactory(fitterFactory) {
  config.setupDefaults();
  parser.getHPOConfig(config);
}
} /* namespace datadriven */
} /* namespace sgpp */
