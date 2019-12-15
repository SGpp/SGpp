// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ClassificationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ClusteringFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#ifdef USE_BOOST_GRAPH
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#endif

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase *UniversalMinerFactory::createFitter(const DataMiningConfigParser &parser) const {
  ModelFittingBase *model = nullptr;

  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);
  if (fType == FitterType::DensityEstimation) {
    FitterConfigurationDensityEstimation config{};
    config.readParams(parser);
#ifdef USE_SCALAPACK
    if (parser.hasParallelConfig()) {
      model = new ModelFittingDensityEstimationOnOffParallel(config);
    } else {
      model = new ModelFittingDensityEstimationOnOff(config);
    }
#else
    model = new ModelFittingDensityEstimationOnOff(config);
#endif /* USE_SCALAPACK */
  } else if (fType == FitterType::RegressionLeastSquares) {
    FitterConfigurationLeastSquares config{};
    config.readParams(parser);
    model = new ModelFittingLeastSquares(config);
  } else if (fType == FitterType::Classification) {
    FitterConfigurationClassification config{};
    config.readParams(parser);
    model = new ModelFittingClassification(config);
  }
  #ifdef USE_BOOST_GRAPH
  else if (fType == FitterType::Clustering) {
    FitterConfigurationClustering config {};
    config.readParams(parser);
    model = new ModelFittingClustering(config);
  }
  #endif
  return model;
}

FitterFactory *UniversalMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);
  FitterFactory *fitfac = nullptr;

  if (fType == FitterType::DensityEstimation) {
    fitfac = new DensityEstimationFitterFactory(parser);
  } else if (fType == FitterType::RegressionLeastSquares) {
    fitfac = new LeastSquaresRegressionFitterFactory(parser);
  } else if (fType == FitterType::Classification) {
    fitfac = new ClassificationFitterFactory(parser);
  }
  #ifdef USE_BOOST_GRAPH
  else if (fType == FitterType::Clustering) {
    fitfac = new ClusteringFitterFactory(parser);
  }
  #endif
  return fitfac;
}

Visualizer *UniversalMinerFactory::createVisualizer(const DataMiningConfigParser &parser) const {
  Visualizer *visualizer = nullptr;

  VisualizerConfiguration config;

  config.readParams(parser);
  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);
  if (fType == FitterType::DensityEstimation) {
    visualizer = new VisualizerDensityEstimation(config);
  } else if (fType == FitterType::RegressionLeastSquares) {
    visualizer = new VisualizerDummy();
  } else if (fType == FitterType::Classification) {
    visualizer = new VisualizerClassification(config);
  }
  #ifdef USE_BOOST_GRAPH
  else if (fType == FitterType::Clustering) {
    visualizer = new VisualizerDummy();
  }
  #endif
  return visualizer;
}
} /* namespace datadriven */
} /* namespace sgpp */
