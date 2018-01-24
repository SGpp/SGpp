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

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource* dataSource, FitterFactory* fitterFactory, Scorer* scorer, HPOScorer* hpoScorer)
    : dataSource(dataSource), fitterFactory(fitterFactory), scorer(scorer), hpoScorer(hpoScorer) {}

void HyperparameterOptimizer::optimizeHyperparameters(){
  // prepare data and scorer
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  double stdDeviation;
  
  // get configured models for n samples (call fitterfactory)
  
  // build matrix
  
  // run samples (parallel)

  // run solver

  //double score = hpoScorer->calculateScore(*fitter, *dataset, &stdDeviation);
}

} /* namespace datadriven */
} /* namespace sgpp */
