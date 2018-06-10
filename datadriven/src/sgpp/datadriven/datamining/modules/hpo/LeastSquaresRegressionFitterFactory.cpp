/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFitterFactory.cpp
 *
 *  Created on:	17.12.2017
 *      Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

LeastSquaresRegressionFitterFactory::LeastSquaresRegressionFitterFactory(DataMiningConfigParser &parser)
    : baseConfig() {

  baseConfig.readParams(parser);

  parser.getHyperparameters(conpar, dispar, catpar, basisFunctions);

  /* dispar["noPoints"] = DiscreteParameter("noPoints",1,4);

   dispar["level"] = DiscreteParameter("level", 1, 4);

     catpar["basisFunction"] = DiscreteParameter("basisFunction",0,1);

     conpar["lambda"] = ContinuousParameter(5, "lambda", -7, 0); //8

     conpar["threshold"] = ContinuousParameter(3, "threshold", -5, -2); //3
 */
}

ModelFittingBase *LeastSquaresRegressionFitterFactory::buildFitter() {

  // build config
  auto *config = new FitterConfigurationLeastSquares(baseConfig);

  //base::GridType basisFunction[] = {base::GridType::Linear,base::GridType::ModLinear};
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
    config->getRefinementConfig().threshold_ = pow(10, conpar["threshold"].getValue());
  }
  if (conpar.count("lambda")) {
    config->getRegularizationConfig().lambda_ = pow(10, conpar["lambda"].getValue());
  }
  return new ModelFittingLeastSquares(*config);
}

// EDIT: exponential parameters

} /* namespace datadriven */
} /* namespace sgpp */
