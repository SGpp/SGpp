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

SparseGridMiner::~SparseGridMiner() {}

void SparseGridMiner::learn() {
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  double stdDeviation;
  double score = scorer->calculateScore(*fitter, *dataset, &stdDeviation);

  std::cout << "Learner finished." << std::endl
            << "Score: " << score << std::endl
            << " Standard Deviation: " << stdDeviation << std::endl;
}

} /* namespace datadriven */
} /* namespace sgpp */
