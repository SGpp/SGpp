/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DensityEstimationFitterFactory.cpp
 *
 *  Created on:	17.12.2017
 *      Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

namespace sgpp {
namespace datadriven {

DensityEstimationFitterFactory::DensityEstimationFitterFactory(DataMiningConfigParser &parser)
    : baseConfig() {

  baseConfig.readParams(parser);

  //EDIT: new hier ohne pointer?

  dispar["level"] = DiscreteParameter("level", 4, 7);

  //catpar["basisFunction"] = DiscreteParameter("basisFunction",0,1);

  conpar["lambda"] = ContinuousParameter(7, "lambda", -10, 0);

}

ModelFittingBase *DensityEstimationFitterFactory::buildFitter() {

  // build config
  auto *config = new FitterConfigurationDensityEstimation(baseConfig);
  //EDIT: make lambda exponential
  base::GridType basisFunction[] = {base::GridType::Linear, base::GridType::ModLinear};

  config->getGridConfig().level_ = dispar["level"].getValue();
  //config->getGridConfig().type_ = basisFunction[catpar["basisFunction"].getValue()];
  config->getRegularizationConfig().lambda_ = pow(10, conpar["lambda"].getValue());

  return new ModelFittingDensityEstimation(*config);
}


} /* namespace datadriven */
} /* namespace sgpp */
