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

#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>


namespace sgpp {
namespace datadriven {



int FitterFactory::buildParity(){
	freeBits = std::vector<ConfigurationBit*>{};
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

int FitterFactory::addConstraint(int idx, int bias){
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
    for(auto& bit : configBits){
      if(!bit->checkConstraints()){
        //EDIT: revert constraint
        for(auto &bit : parityrow[idx]){
          bit->removeLastConstraint();
          std::cout<<"Adding bit from constraint:"<<idx<<std::endl;
        }
      }
    }
		return freeBits.size();
}

/*EDIT: check old runs in new space: set free bits, evaluate, then check constraints
-error handling: prevent conflicting constraints
 -handle constraints with zero bias
*/

void FitterFactory::setHarmonica(int configID, int row, DataMatrix &paritymatrix)  {
  // fix ConfigurationBits according to constraints
  for(auto& bit : configBits){
    bit->reset();
  }

  for(auto& bit : freeBits){
    bit->setBit(&configID);
  }
  // remove, do this through parameters
	// for(auto& bit : configBits){
  //  bit.evaluate(&configID);
	// }


  for(int i=0;i<parityrow.size();i++){
	  int tmp = 1;
	  for(auto& bit : parityrow[i]){
		  tmp = tmp * bit->evaluate();
	  }
	  //std::cout<<"Bitinparity:"<<tmp<<std::endl;
	  paritymatrix.set(row, i, tmp);
  }

  if(configID != 0){
    std::cout<<"Error: configID not fully used:"<<configID<<std::endl;
  }

  for(auto pair: dispar){
    pair.second->setHarmonica();
  }

  for(auto pair: conpar){
    pair.second->setHarmonica();
  }
  // std::cout<<"Run mark 2.1"<<std::endl;
}

void FitterFactory::getBOspace(int* nCont, std::vector<int>& nOptions){
  *nCont = conpar.size();
  for(auto pair: dispar){
    nOptions.push_back(pair.second->getNOptions());
  }
} //EDIT: add categorical parameters

void FitterFactory::setBO(base::DataVector& cont, std::vector<int>& disc){
  int i = 0;
  for(auto pair: dispar){
    pair.second->setBO(disc[i]);
    i++;
  }
  i = 0;
  for(auto pair: conpar){
    pair.second->setBO(cont[i]);
    i++;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
