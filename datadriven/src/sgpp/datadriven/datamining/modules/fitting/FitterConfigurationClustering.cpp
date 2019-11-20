// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClustering.hpp>

namespace sgpp {
namespace datadriven {

  FitterConfiguration* FitterConfigurationClustering::clone() const {
    return new FitterConfigurationClustering(*this);
  }

  void FitterConfigurationClustering::setupDefaults() {
    FitterConfiguration::setupDefaults();
  }

  void FitterConfigurationClustering::readParams(const DataMiningConfigParser &parser) {
    setupDefaults();
    FitterConfigurationDensityEstimation::readParams(parser);
    parser.getFitterClusteringConfig(clusteringConfig, clusteringConfig);
  }
}  // namespace datadriven
}  // namespace sgpp
