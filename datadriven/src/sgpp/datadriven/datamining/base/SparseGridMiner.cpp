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
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

SparseGridMiner::SparseGridMiner(ModelFittingBase* fitter, Scorer* scorer)
  : fitter(fitter), scorer(scorer) {}

double SparseGridMiner::test(Dataset& testDataset) {
  return scorer->test(*fitter, testDataset);
}

ModelFittingBase *SparseGridMiner::getModel() {
  return &(*fitter);
}

void SparseGridMiner::setModel(ModelFittingBase *model) {
  fitter.reset(model);
}
} /* namespace datadriven */
} /* namespace sgpp */
