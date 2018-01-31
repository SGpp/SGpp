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

#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

LeastSquaresRegressionFitterFactory::LeastSquaresRegressionFitterFactory(DataMiningConfigParser& parser)
  :configBits(), baseConfig(), parityrow(){
  // build ConfigurationBits (constructor)
  // for(int i=0;i<12;i++){
  //  configBits.push_back(*(new ConfigurationBit()));
  // }
	baseConfig.readParams(parser);

	//EDIT: new hier ohne pointer?
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  conpar["lambda"] = (new ContinuousParameter{-2, 1});
  conpar["lambda"]->makeConfigBits(4, configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  conpar["threshold"] = (new ContinuousParameter{0.0005, 0.002});
  conpar["threshold"]->makeConfigBits(3, configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  dispar["noPoints"] = (new DiscreteParameter{1,4});
  dispar["noPoints"]->makeConfigBits(configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

  dispar["level"] = (new DiscreteParameter{1,4});
  dispar["level"]->makeConfigBits(configBits);
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;

	dispar["basisFunction"] = (new DiscreteParameter{0,1});
	dispar["basisFunction"]->makeConfigBits(configBits);

}

int LeastSquaresRegressionFitterFactory::buildParity(){
	std::vector<ConfigurationBit*> freeBits{};
	for(auto& bit : configBits){ //reference to prevent unique pointer copying
	    bit->reset();
	  }
	for(auto& bit : configBits){
		bit->fixFreeBits(freeBits);
	}
	int ncols = freeBits.size();
	std::cout<<"nBits: "<<ncols<<std::endl;
	std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;
	parityrow.clear();
	parityrow.reserve((ncols*ncols+5)*ncols/6 +1);
	for(int i = 0; i<(ncols*ncols+5)*ncols/6 +1;i++){
		parityrow.push_back(*(new std::list<ConfigurationBit*>()));
	}
	int cnt = ncols;
	int cnt2 =(ncols +1)*ncols/2;
	  std::cout<<"Run mark 1.1"<<std::endl;
	for(int i = 0; i < ncols; i++){
		//parityrow[i] = *new std::list<ConfigurationBit*>();
		  //std::cout<<"Run mark 1.1:"<<i<<std::endl;
		parityrow[i].push_back(freeBits[i]);
		for(int k = i+1; k < ncols; k++){
		//	parityrow[cnt] = *new std::list<ConfigurationBit*>();
			parityrow[cnt].push_back(freeBits[i]);
			parityrow[cnt].push_back(freeBits[k]);
			if(cnt==77){
				  std::cout<<"Constraint 297: "<<i<<","<<k<<std::endl;
			}
	        cnt++;
	        for(int m = k+1; m < ncols; m++){
		//		parityrow[cnt2] = *new std::list<ConfigurationBit*>();
				parityrow[cnt2].push_back(freeBits[i]);
				parityrow[cnt2].push_back(freeBits[k]);
				parityrow[cnt2].push_back(freeBits[m]);
				if(cnt2==297){
					  std::cout<<"Constraint 297: "<<i<<","<<k<<","<<m<<std::endl;
				}
				cnt2++;
			}
		}
	}
	  std::cout<<"Run mark 1.2"<<std::endl;
	return ncols;
}

int LeastSquaresRegressionFitterFactory::addConstraint(int idx, int bias){
	ConfigurationRestriction* hey = new ConfigurationRestriction(parityrow[idx], bias); //EDIT: store and delete
	for(auto &bit : parityrow[idx]){
		bit->addConstraint(hey);
		std::cout<<"Adding bit from constraint:"<<idx<<std::endl;
	}
	std::vector<ConfigurationBit*> freeBits{};
		for(auto& bit : configBits){ //reference to prevent unique pointer copying
		    bit->reset();
		  }
		for(auto& bit : configBits){
			bit->fixFreeBits(freeBits);
		}
		return freeBits.size();
}

//check old runs in new space: set free bits, evaluate, then check constraints

void LeastSquaresRegressionFitterFactory::printConfig(int configID){
	  for(auto& bit : configBits){
	    bit->reset();
	  }
	  std::string basisFunction[] = {"Linear","ModLinear"};
	  std::cout<<"Level: "<<dispar["level"]->getValue(&configID)
			   <<", Basis: "<<basisFunction[dispar["basisFunction"]->getValue(&configID)]
			   <<", noPoints: "<< dispar["noPoints"]->getValue(&configID)
			   <<", Threshold: "<<conpar["threshold"]->getValue(&configID)
			   <<", Lambda: "<<conpar["lambda"]->getValue(&configID)<<", "<<pow(10,conpar["lambda"]->getValue(&configID))
			   <<std::endl;
}

ModelFittingBase* LeastSquaresRegressionFitterFactory::buildFitter(int configID, int row, DataMatrix &paritymatrix)  {
  // fix ConfigurationBits according to constraints
  for(auto& bit : configBits){
    bit->reset();
  }
  // remove, do this through parameters
  // for(auto bit : configBits){
  //  bit.evaluate(&configID);
  // }

  // build config
  FitterConfigurationLeastSquares* config = new FitterConfigurationLeastSquares(baseConfig);
  //EDIT: make lambda exponential
  base::GridType basisFunction[] = {base::GridType::Linear,base::GridType::ModLinear};
  config->setHyperParameters(
		  dispar["level"]->getValue(&configID),
		  basisFunction[dispar["basisFunction"]->getValue(&configID)],
		  dispar["noPoints"]->getValue(&configID),
		  conpar["threshold"]->getValue(&configID),
		  pow(10,conpar["lambda"]->getValue(&configID))
		  );
 /* std::cout<<"Error: configID not fully used:"<<configID<<std::endl;

  dispar["level"]->getValue(&configID);
  std::cout<<"Error: configID not fully used:"<<configID<<std::endl;

  dispar["noPoints"]->getValue(&configID);
  std::cout<<"Error: configID not fully used:"<<configID<<std::endl;

  conpar["threshold"]->getValue(&configID);
  std::cout<<"Error: configID not fully used:"<<configID<<std::endl;

		conpar["lambda"]->getValue(&configID);
		  std::cout<<"last: configID not fully used:"<<configID<<std::endl;
*/
  // return model and ConfigurationBits in vector/matrix
  if(configID != 0){
	  std::cout<<"Error: configID not fully used:"<<configID<<std::endl;
  }
  for(int i=0;i<parityrow.size();i++){
	  int tmp = 1;
	  for(auto& bit : parityrow[i]){
		  tmp = tmp * bit->evaluate(&configID);
	  }
	  //std::cout<<"Bitinparity:"<<tmp<<std::endl;
	  paritymatrix.set(row, i, tmp);
  }
  // std::cout<<"Run mark 2.1"<<std::endl;

  return new ModelFittingLeastSquares(*config);
}
FitterConfiguration* LeastSquaresRegressionFitterFactory::buildConfig() const {
  FitterConfigurationLeastSquares* config = new FitterConfigurationLeastSquares(baseConfig);
  // DataMiningConfigParser parser("TrainWithTestingExample.json");
  // config->readParams(parser);

  return config;
}
} /* namespace datadriven */
} /* namespace sgpp */
