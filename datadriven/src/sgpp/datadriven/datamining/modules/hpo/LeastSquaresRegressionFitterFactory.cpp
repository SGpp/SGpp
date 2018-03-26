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

LeastSquaresRegressionFitterFactory::LeastSquaresRegressionFitterFactory(DataMiningConfigParser& parser)
  :baseConfig(){

	baseConfig.readParams(parser);

  dispar["noPoints"] = DiscreteParameter("noPoints",1,4);

  dispar["level"] = DiscreteParameter("level", 1, 4);

	catpar["basisFunction"] = DiscreteParameter("basisFunction",0,1);

	conpar["lambda"] = ContinuousParameter(5, "lambda", -7, 0); //8

	conpar["threshold"] = ContinuousParameter(3, "threshold", -5, -2); //3

}




ModelFittingBase* LeastSquaresRegressionFitterFactory::buildFitter()  {

  // build config
  auto* config = new FitterConfigurationLeastSquares(baseConfig);

  base::GridType basisFunction[] = {base::GridType::Linear,base::GridType::ModLinear};
  config->getGridConfig().level_ = dispar["level"].getValue();
  config->getGridConfig().type_ = basisFunction[catpar["basisFunction"].getValue()];
  config->getRefinementConfig().noPoints_ = static_cast<size_t>(dispar["noPoints"].getValue());
  config->getRefinementConfig().threshold_ = pow(10,conpar["threshold"].getValue());
  config->getRegularizationConfig().lambda_ = pow(10,conpar["lambda"].getValue());

  return new ModelFittingLeastSquares(*config);
}

std::string LeastSquaresRegressionFitterFactory::printConfig(){

	std::string basisFunction[] = {"Linear","ModLinear"};
	std::cout<<"Level: "<< dispar["level"].getValue()
					 <<", Basis: "<<basisFunction[catpar["basisFunction"].getValue()]
					 <<", noPoints: "<< dispar["noPoints"].getValue()
					 <<", Threshold: "<< pow(10,conpar["threshold"].getValue())
					 <<", Lambda: "<< pow(10,conpar["lambda"].getValue())
					 <<std::endl;
  return std::to_string(dispar["level"].getValue())
                    + ", " + basisFunction[catpar["basisFunction"].getValue()]
                    + ", " + std::to_string(dispar["noPoints"].getValue())
                    + ", " + std::to_string(pow(10,conpar["threshold"].getValue()))
                    + ", " + std::to_string(pow(10,conpar["lambda"].getValue()));
                   // + std::endl;
}

} /* namespace datadriven */
} /* namespace sgpp */
