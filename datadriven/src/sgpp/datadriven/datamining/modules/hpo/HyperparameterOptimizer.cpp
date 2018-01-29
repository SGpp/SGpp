/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HyperparameterOptimizer.cpp
 *
 * Created on: Jan 22, 2018
 *     Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/OperationMultipleEvalMatrix.hpp>
#include <sgpp/base/Grid/type/LinearGrid.hpp>


#include <iostream>
#include <random>
#include <cmath>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource* dataSource, FitterFactory* fitterFactory, Scorer* scorer, HPOScorer* hpoScorer)
    : dataSource(dataSource), fitterFactory(fitterFactory), scorer(scorer), hpoScorer(hpoScorer) {}

void HyperparameterOptimizer::optimizeHyperparameters(){
  // prepare data and scorer
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());

  std::cout<<"Run mark 1"<<std::endl;

  
  // get configured models for n samples (call fitterfactory)
  int nBits = fitterFactory->buildParity();
  int n = 50;
  
  //build matrix
  // EDIT: add bias term?
  DataMatrix paritymatrix(n, (nBits*nBits+5)*nBits/6 +1);

  //make random inputs
  int seed = 42;
  std::mt19937 generator = std::mt19937(seed);
  std::uniform_int_distribution<int> distribution(0,std::pow(2,nBits));
  std::cout<<"MaxConfigs: "<<std::pow(2,nBits)<<std::endl;
  std::vector<int> configIDs(n);
  for(int i=0;i<n;i++){
	  configIDs[i] = distribution(generator);
	  bool bUnchecked = true;
	  while(bUnchecked){
		  bUnchecked = false;
		  for(int k = 0; k<i; k++){
			  if(configIDs[i]==configIDs[k]){
			  	  configIDs[i] = distribution(generator);
			  	  bUnchecked = true;
			  	  //std::cout<<"Double on Config: "<<k<<";"<<configIDs[k]<<std::endl;
		  	  }
	  	  }
	  }
  }
  std::cout<<"Run mark 2"<<std::endl;

  std::vector<ModelFittingBase*> fitters(n);
  for(int i=0;i<n;i++){
	  fitters[i] = fitterFactory->buildFitter(configIDs[i], i, paritymatrix);
  }
  
  std::cout<<"Run mark 2.5"<<std::endl;

  // run samples (parallel)
  DataVector scores(n);
  double stdDeviation;
  for(int i=0;i<n;i++){
	  scores[i] = hpoScorer->calculateScore(*(fitters[i]), *dataset, &stdDeviation);
	  std::cout<<"Score "<<i<<":"<<scores[i]<<std::endl;
  }

  std::cout<<"Run mark 3"<<std::endl;

  std::cout<< "size: "<<paritymatrix.getNrows()<<" x "<<paritymatrix.getNcols()<<", "<<paritymatrix.get(10, 100)<<", "<<paritymatrix.get(10, 10)<<std::endl;


  // run solver
  solver::LassoFunction g{10}; //EDIT: lambda
  solver::Fista<solver::LassoFunction> fista{g};
  DataVector alpha = DataVector{paritymatrix.getNcols()};
  OperationMultipleEvalMatrix opMultEval{*new base::LinearGrid(0), paritymatrix}; //EDIT: grid löschen
  fista.solve(opMultEval, alpha, scores, 100, DEFAULT_RES_THRESHOLD);

  std::vector<int> idx(alpha.size());
  for(int i=0; i<idx.size(); i++){
	  idx[i] = i;
	  std::cout<<"Alpha: "<<i<<":"<<alpha[i]<<std::endl;
  }

   // sort indices based on comparing values in alpha
   sort(idx.begin(), idx.end(),
         [&alpha](int i1, int i2) {return fabs(alpha[i1]) > fabs(alpha[i2]);});

   int shrink = 1; //EDIT: parameter

   for(int i = 0; i < shrink; i++){
	   fitterFactory->addConstraint(idx[i], lround(alpha[idx[i]]/fabs(alpha[idx[i]]))); //(alpha[idx[i]]>0)-(alpha[idx[i]]<0)
	   std::cout<<"Constraint number:"<<idx[i]<<","<<lround(alpha[idx[i]]/fabs(alpha[idx[i]]))<<","<<alpha[idx[i]]<<std::endl;
   }

   std::cout<<"Run mark 4"<<std::endl;




    nBits = fitterFactory->buildParity();
      n = 64;

     //build matrix
     // EDIT: add bias term?
     paritymatrix.resize(n, (nBits*nBits+5)*nBits/6 +1);

     //make random inputs
     //int seed = 42;
     //std::mt19937 generator = std::mt19937(seed);
     std::uniform_int_distribution<int> distribution2(0,std::pow(2,nBits));
     std::cout<<"MaxConfigs: "<<std::pow(2,nBits)<<std::endl;
     //std::vector<int> configIDs2(n);
     for(int i=0;i<n;i++){
   	  configIDs[i] = distribution2(generator);
   	  bool bUnchecked = true;
   	  while(bUnchecked){
   		  bUnchecked = false;
   		  for(int k = 0; k<i; k++){
   			  if(configIDs[i]==configIDs[k]){
   			  	  configIDs[i] = distribution2(generator);
   			  	  bUnchecked = true;
   			  	  //std::cout<<"Double on Config: "<<k<<";"<<configIDs[k]<<std::endl;
   		  	  }
   	  	  }
   	  }
     }
     std::cout<<"Run mark 2"<<std::endl;

     //std::vector<ModelFittingBase*> fitters(n);
     for(int i=0;i<n;i++){
   	  fitters[i] = fitterFactory->buildFitter(configIDs[i], i, paritymatrix);
     }

     std::cout<<"Run mark 2.5"<<std::endl;

     // run samples (parallel)
     //DataVector scores(n);
     //double stdDeviation;
     for(int i=0;i<n;i++){
   	  scores[i] = hpoScorer->calculateScore(*(fitters[i]), *dataset, &stdDeviation);
   	  std::cout<<"Score "<<i<<":"<<scores[i]<<std::endl;
     }





  //double score = hpoScorer->calculateScore(*fitter, *dataset, &stdDeviation);
}

} /* namespace datadriven */
} /* namespace sgpp */
