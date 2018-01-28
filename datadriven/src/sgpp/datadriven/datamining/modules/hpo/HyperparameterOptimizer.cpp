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
#include <sgpp/datadriven/tools/CSVTools.hpp>


#include <iostream>
#include <random>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource* dataSource, FitterFactory* fitterFactory, Scorer* scorer, HPOScorer* hpoScorer)
    : dataSource(dataSource), fitterFactory(fitterFactory), scorer(scorer), hpoScorer(hpoScorer) {}

void HyperparameterOptimizer::optimizeHyperparameters(){
  // prepare data and scorer
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  double stdDeviation;
  
  // get configured models for n samples (call fitterfactory)
  int nBits = fitterFactory.buildParity();
  int n = 300;
  
  //build matrix
  // EDIT: add bias term?
  DataMatrix paritymatrix(n, (nBits*nBits+5)*nBits/6 +1);

  //make random inputs
  int seed = 42;
  std::mt19937 generator = std::mt19937(seed);
  std::uniform_int_distribution<int> distribution(0,2**nBits-1);
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
		  	  }
	  	  }
	  }
  }
  std::vector<ModelFittingBase*> fitters{};
  for(int i=0;i<n;i++){
	  ModelFittingBase* fitters[i] = fitterFactory.buildFitter(configIDs[i], n, paritymatrix);
  }
  

  // run samples (parallel)
  DataVector scores(n);
  for(int i=0;i<n;i++){
	  scores[i] = hpoScorer->calculateScore(*fitters[i], *dataset, &stdDeviation);
  }


  // run solver
  solver::LassoFunction g{100}; //EDIT: lambda
  solver::Fista<solver::LassoFunction> fista{g};
  alpha = DataVector{paritymatrix.getNcols()};
  OperationMultipleEvalMatrix opMultEval{*grid, paritymatrix}; //EDIT: grid raustun
  fista.solve(opMultEval, alpha, scores, 100, DEFAULT_RES_THRESHOLD);

  std::vector<int> idx(alpha.size());
  for(int i=0; i<idx.size(); i++){
	  idx[i] = i;
  }

   // sort indices based on comparing values in alpha
   sort(idx.begin(), idx.end(),
         [&v](int i1, int i2) {return abs(alpha[i1]) > abs(alpha[i2]);});

   int shrink = 5; //EDIT: parameter

   for(int i = 0; i < shrink; i++){
	   fitterFactory->addConstraint(idx[i], lround(alpha[idx[i]]/abs(alpha[idx[i]])));
   }

  //double score = hpoScorer->calculateScore(*fitter, *dataset, &stdDeviation);
}

} /* namespace datadriven */
} /* namespace sgpp */
