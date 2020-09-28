// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/RelativeDensityRatioEstimationMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplittingTwoDatasets.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingRelativeDensityRatioEstimation.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

SparseGridMiner* RelativeDensityRatioEstimationMinerFactory::buildMiner(
    const std::string& path) const {
  DataMiningConfigParser parser(path);
  if (parser.hasFitterConfigCrossValidation()) {
    // TODO(fuchsgdk): implement the cv stuff
    return new SparseGridMinerCrossValidation(createDataSourceCrossValidation(parser),
                                              createFitter(parser), createScorer(parser),
                                              createVisualizer(parser));
  } else {
    return new SparseGridMinerSplittingTwoDatasets(createDataSourceSplittingTwoDatasets(parser),
                                                   createFitter(parser), createScorer(parser),
                                                   createVisualizer(parser));
  }
}

std::vector<DataSourceSplitting*>
RelativeDensityRatioEstimationMinerFactory::createDataSourceSplittingTwoDatasets(
    const DataMiningConfigParser& parser) const {
  std::vector<DataSourceConfig> configs(2);

  bool hasSource = parser.getMultiDataSourceConfig(configs, configs);

  std::vector<DataSourceSplitting*> dataSources(2);

  for (size_t i = 0; i < dataSources.size(); ++i)
    if (hasSource && configs[i].filePath_.compare("") != 0) {
      DataSourceBuilder builder;
      dataSources[i] = builder.splittingFromConfig(configs[i]);
    } else {
      throw base::data_exception("No file name provided for datasource.");
    }

  return dataSources;
}

/*
 std::vector<DataSourceCrossValidation*>
 RelativeDensityRatioEstimationMinerFactory::createDataSourceCrossValidationTwoDatasets(
 const DataMiningConfigParser& parser) const {
 DataSourceConfig config { };

 CrossvalidationConfiguration crossValidationconfig { };
 parser.getFitterCrossvalidationConfig(crossValidationconfig,
 crossValidationconfig);

 bool hasSource = parser.getDataSourceConfig(config, config);

 if (hasSource && config.filePath.compare("") != 0) {
 DataSourceBuilder builder;
 return builder.crossValidationFromConfig(config, crossValidationconfig);
 } else {
 throw base::data_exception("No file name provided for datasource.");
 }
 }
 */

ModelFittingBase* RelativeDensityRatioEstimationMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  return new ModelFittingRelativeDensityRatioEstimation(config);
}

/*
 FitterFactory* RelativeDensityRatioEstimationMinerFactory::createFitterFactory(
 const DataMiningConfigParser& parser) const {
 return new DensityRatioEstimationMinerFactory(parser);
 }
 */

Visualizer* RelativeDensityRatioEstimationMinerFactory::createVisualizer(
    const DataMiningConfigParser& parser) const {
  VisualizerConfiguration config;
  config.readParams(parser);

  // TODO(spc90): implement visualization for this model
  return new VisualizerDummy();
}

} /* namespace datadriven */
} /* namespace sgpp */
