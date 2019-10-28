// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/ClusteringFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>

namespace sgpp {
namespace datadriven {

  ClusteringFitterFactory::ClusteringFitterFactory
      (const DataMiningConfigParser &parser) : baseConfig() {
    baseConfig.readParams(parser);

    parser.getHyperparameters(conpar, dispar, catpar, basisFunctions);
  }

  ModelFittingBase *ClusteringFitterFactory::buildFitter() {
    // build config
    auto *config = new FitterConfigurationClustering(baseConfig);

    if (dispar.count("level")) {
      config->getGridConfig().level_ = dispar["level"].getValue();
    }
    if (catpar.count("basisFunction")) {
      config->getGridConfig().type_ = basisFunctions[catpar["basisFunction"].getValue()];
    }
    if (dispar.count("noPoints")) {
      config->getRefinementConfig().noPoints_ = static_cast<size_t>(dispar["noPoints"].getValue());
    }
    if (conpar.count("threshold")) {
      config->getRefinementConfig().threshold_ = conpar["threshold"].getValue();
    }
    if (conpar.count("lambda")) {
      config->getRegularizationConfig().lambda_ = conpar["lambda"].getValue();
    }

    return new ModelFittingClustering(*config);
  }
}  // namespace datadriven
}  // namespace sgpp
