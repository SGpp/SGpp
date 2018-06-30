/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.cpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

SparseGridMiner::SparseGridMiner(DataSource* dataSource, ModelFittingBase* fitter, Scorer* scorer)
    : dataSource(dataSource), fitter(fitter), scorer(scorer) {}

void SparseGridMiner::learn() {
  double totalScore = 0.0, totalStdDeviation = 0.0;

  // Process the dataset iteratively
  size_t iteration = 0;
  while (true) {
    std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
    size_t numInstances = dataset->getNumberInstances();
    if (numInstances == 0) {
      // The source does not provide any more samples
      break;
    }
    std::cout <<  "###############" << "Dataset iteration " << (iteration++) << std::endl <<
        "Number of instances: " << numInstances << std::endl;
    double stdDeviation;
    double score = scorer->calculateScore(*fitter, *dataset, &stdDeviation);
    std::cout << "Iteration finished." << std::endl
              << "###############" << std::endl
              << "Score: " << score << std::endl
              << "Standard Deviation: " << stdDeviation << std::endl
              << "###############" << std::endl;
    totalScore += score;
    totalStdDeviation += stdDeviation;
  }
  totalScore /= static_cast<double>(iteration);
  totalStdDeviation /= static_cast<double>(iteration);
  std::cout << "###############" << std::endl << "Learner finished." << std::endl <<
      "Mean score: " << totalScore << std::endl << "Mean standard deviation: " <<
      totalStdDeviation << std::endl;
}

ModelFittingBase *SparseGridMiner::getModel() {
  return &(*fitter);
}

} /* namespace datadriven */
} /* namespace sgpp */
