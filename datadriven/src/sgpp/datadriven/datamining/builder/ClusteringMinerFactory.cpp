// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/ClusteringMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ClusteringFitterFactory.hpp>

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClustering.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

    ModelFittingBase* ClusteringMinerFactory::createFitter(
        const DataMiningConfigParser& parser) const {
      FitterConfigurationClustering config{};
      config.readParams(parser);
      return new ModelFittingClustering(config);
    }

    FitterFactory *ClusteringMinerFactory::createFitterFactory(
        const DataMiningConfigParser &parser) const {
      return new ClusteringFitterFactory(parser);
    }

    Visualizer* ClusteringMinerFactory::createVisualizer(const DataMiningConfigParser& parser)
    const {
      VisualizerConfiguration config;

      config.readParams(parser);

      return new VisualizerClustering(config);
    }

} /* namespace datadriven */
} /* namespace sgpp */