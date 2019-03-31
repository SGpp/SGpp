/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * UniversalMinerFactory.cpp
 *
 * Created on: Mar 12, 2018
 *     Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include <string>

namespace sgpp {
namespace datadriven {


ModelFittingBase *UniversalMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  ModelFittingBase *model;

  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);
  if (fType == FitterType::DensityEstimation) {
    FitterConfigurationDensityEstimation config{};
    config.readParams(parser);
    if (!config.getCombi()) {
      model = new ModelFittingDensityEstimationOnOff(config);
    } else {
       return new PDFCombigrid(config);
    }
  } else {
    FitterConfigurationLeastSquares config{};
    config.readParams(parser);
    model = new ModelFittingLeastSquares(config);
  }
  return model;
}

FitterFactory *UniversalMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);
  FitterFactory* fitfac;

  if (fType == FitterType::DensityEstimation) {
    fitfac = new DensityEstimationFitterFactory(parser);
  } else {
    fitfac = new LeastSquaresRegressionFitterFactory(parser);
  }
  return fitfac;
}


} /* namespace datadriven */
} /* namespace sgpp */
