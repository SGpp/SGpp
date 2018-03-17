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
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

DensityEstimationFitterFactory::DensityEstimationFitterFactory(DataMiningConfigParser& parser)
  :baseConfig(){

	baseConfig.readParams(parser);

	//EDIT: new hier ohne pointer?

  dispar["level"] = (new DiscreteParameter{1,6});
  dispar["level"]->makeConfigBits(configBits);

  //EDIT: catpar here
  dispar["basisFunction"] = (new DiscreteParameter{0,1});
  dispar["basisFunction"]->makeConfigBits(configBits);

	conpar["lambda"] = (new ExponentialParameter{-8, 0});
	conpar["lambda"]->makeConfigBits(4, configBits);

}



ModelFittingBase* DensityEstimationFitterFactory::buildFitter()  {

  // build config
  FitterConfigurationLeastSquares* config = new FitterConfigurationLeastSquares(baseConfig);
  //EDIT: make lambda exponential
  base::GridType basisFunction[] = {base::GridType::Linear,base::GridType::ModLinear};

  config->getGridConfig().level_ = dispar["level"]->getValue();
  config->getGridConfig().type_ = basisFunction[dispar["basisFunction"]->getValue()];
  config->getRegularizationConfig().lambda_ = conpar["lambda"]->getValue();

  return new ModelFittingLeastSquares(*config);
}

void DensityEstimationFitterFactory::printConfig(){

	std::string basisFunction[] = {"Linear","ModLinear"};
	std::cout<<"Level: "<< dispar["level"]->getValue()
					 <<", Basis: "<<basisFunction[dispar["basisFunction"]->getValue()]
					 <<", Lambda: "<<conpar["lambda"]->getValue()
					 <<std::endl;
}

} /* namespace datadriven */
} /* namespace sgpp */
