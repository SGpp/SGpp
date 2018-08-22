// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/BatchLearner.hpp"
#include "sgpp/datadriven/configuration/BatchConfiguration.hpp"
#include "sgpp/datadriven/tools/Dataset.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"

/**
 * This programm demonstrates the usage of the BatchLearner class. After the
 * parameters are set, the method trainBatch() is called until the end of the
 * file has been reached.
 */

int main(int argc, char** args) {
  // set variables
  sgpp::base::BatchConfiguration batchConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::base::AdaptivityConfiguration adaptConfig;
  sgpp::base::RegularGridConfiguration gridConfig;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;  // not used by BatchLearner
  adaptConfig.noPoints_ = 2;
  adaptConfig.numRefinements_ = 1;  // not used by BatchLearner
  adaptConfig.percent_ = 10.0;      // not used by BatchLearner
  adaptConfig.threshold_ = 0.001;

  // Set solver during refinement
  solverConfig.eps_ = 0;
  solverConfig.maxIterations_ = 20;
  solverConfig.threshold_ = -1.0;  // not used by BatchLearner
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;

  // Set parameters for the batchLearner
  batchConfig.filename_ = "../tests/data/friedman_4d_2000.arff";
  batchConfig.batchsize_ = 500;
  batchConfig.samples_ = 500;
  batchConfig.seed_ = 42;
  batchConfig.wMode_ = 5;
  batchConfig.wArgument_ = 1.0;
  batchConfig.refineEvery_ = 0;
  batchConfig.verbose_ = true;
  batchConfig.stack_ = 0;
  batchConfig.testsize_ = 200;
  batchConfig.lambda_ = 0.0001f;

  // set up the grid config
  gridConfig.level_ = 4;

  // init the learner
  sgpp::datadriven::BatchLearner learner(batchConfig, gridConfig, solverConfig,
                                         adaptConfig);

  int i = 0;

  while (!learner.getIsFinished()) {  // while there are still batches to read
    std::cout << "called trainBatch() #" << i << std::endl;
    int64_t start_l = clock();
    learner.trainBatch();
    int64_t stop_l = clock();
    std::cout << "trainBatch() took: "
              << double(stop_l - start_l) / double(CLOCKS_PER_SEC) * 1000
              << std::endl;

    // If testsize > 0, the learner will test on testsize items following the
    // batch; in the next batch these items will be included and the batch is
    // filled with new items until batchsize is reached. This will be done after
    // every batch.
    // AC = accuracy of the current batch
    // AG = accuracy over all predictions
    std::cout << "AC: " << i << ", " << learner.getAccCurrent() << "\nAG: " << i
              << ", " << learner.getAccGlobal() << std::endl;
    i++;
  }

  // don't forget to close the file (only needed if the file has not been read
  // to completion yet, otherwise the reader will close the file automatically)
  learner.closeFile();

  return EXIT_SUCCESS;
}
