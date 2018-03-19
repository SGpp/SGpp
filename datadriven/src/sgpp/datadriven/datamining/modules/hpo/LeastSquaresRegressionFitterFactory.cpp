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
  // build ConfigurationBits (constructor)
  // for(int i=0;i<12;i++){
  //  configBits.push_back(*(new ConfigurationBit()));
  // }
	baseConfig.readParams(parser);

	//EDIT: new hier ohne pointer?
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  dispar["noPoints"] = (new DiscreteParameter{1,4});
  dispar["noPoints"]->makeConfigBits(configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  dispar["level"] = (new DiscreteParameter{1,4});
  dispar["level"]->makeConfigBits(configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

	dispar["basisFunction"] = (new DiscreteParameter{0,1});
	dispar["basisFunction"]->makeConfigBits(configBits);

	conpar["lambda"] = (new ExponentialParameter{-8, 0});
	conpar["lambda"]->makeConfigBits(6, configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

	conpar["threshold"] = (new ExponentialParameter{-5, -2});
	conpar["threshold"]->makeConfigBits(3, configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;



}




ModelFittingBase* LeastSquaresRegressionFitterFactory::buildFitter()  {

  // build config
  FitterConfigurationLeastSquares* config = new FitterConfigurationLeastSquares(baseConfig);
  //EDIT: make lambda exponential
  base::GridType basisFunction[] = {base::GridType::Linear,base::GridType::ModLinear};
  config->getGridConfig().level_ = dispar["level"]->getValue();
  config->getGridConfig().type_ = basisFunction[dispar["basisFunction"]->getValue()];
  config->getRefinementConfig().noPoints_ = dispar["noPoints"]->getValue();
  config->getRefinementConfig().threshold_ = conpar["threshold"]->getValue();
  config->getRegularizationConfig().lambda_ = conpar["lambda"]->getValue();


  return new ModelFittingLeastSquares(*config);
}

void LeastSquaresRegressionFitterFactory::printConfig(){

	std::string basisFunction[] = {"Linear","ModLinear"};
	std::cout<<"Level: "<< dispar["level"]->getValue()
					 <<", Basis: "<<basisFunction[dispar["basisFunction"]->getValue()]
					 <<", noPoints: "<< dispar["noPoints"]->getValue()
					 <<", Threshold: "<<conpar["threshold"]->getValue()<<", "<<pow(10,conpar["threshold"]->getValue())
					 <<", Lambda: "<<conpar["lambda"]->getValue()<<", "<<pow(10,conpar["lambda"]->getValue())
					 <<std::endl;
}

} /* namespace datadriven */
} /* namespace sgpp */
