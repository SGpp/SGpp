/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMinerCrossValidation.cpp
 *
 *  Created on: Jul 26, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

SparseGridMinerCrossValidation::SparseGridMinerCrossValidation(
    DataSourceCrossValidation* dataSource,
    ModelFittingBase* fitter, Scorer* scorer) : SparseGridMiner(fitter, scorer),
        dataSource{dataSource} {}

void SparseGridMinerCrossValidation::learn() {
  // todo(fuchsgdk): see below
  // Setup refinement monitor
  // RefinementMonitorFactory monitorFactory;
  // RefinementMonitor *monitor = monitorFactory.createRefinementMonitor(
  //    fitter->getFitterConfiguration().getRefinementConfig());

  const CrossvalidationConfiguration& crossValidationConfig = dataSource->getCrossValidationConfig();

  double sumScoreVal = 0.0;
  for (size_t fold = 0; fold < crossValidationConfig.kfold_; fold++) {

    // todo(fuchsgdk):
    // This is the kind of cv implemented by Lettrich in the scorer class and it was
    // merely moved to fit into the data source. Conceptual changes might be done in order to
    // really support batch based learning with cv and not only regression.
    // What should be done is reimplementing the data source such that it provides batches

    std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
    std::unique_ptr<Dataset> validation(dataSource->getValidationData());
    size_t numInstances = dataset->getNumberInstances();
    std::cout << "###############" << "Fold #" << fold << std::endl << "Size training data: "
        << numInstances << std::endl;

    // Call the fitting method as it resets the model state
    fitter->fit(*dataset);


  }
}
} /* namespace datadriven */
} /* namespace sgpp */








