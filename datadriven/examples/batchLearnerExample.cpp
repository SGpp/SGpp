// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>
#include <stdio.h>
#include <string>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/BatchLearner.hpp>
#include <sgpp/datadriven/application/BatchConfiguration.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

/**
 * This programm demonstrates the usage of the BatchLearner class. After the parameters are set, the method trainBatch() is called until the end of the file has been reached.
 */

using namespace SGPP::base;
using namespace SGPP::datadriven;
using namespace std;


int main (int argc, char** args) {
  //set variables
  SGPP::base::BatchConfiguration batchConfig;
  SGPP::solver::SLESolverConfiguration solverConfig;
  SGPP::base::AdpativityConfiguration adaptConfig;
  SGPP::base::RegularGridConfiguration gridConfig;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;//not used by BatchLearner
  adaptConfig.noPoints_ = 2;
  adaptConfig.numRefinements_ = 1;//not used by BatchLearner
  adaptConfig.percent_ = 10.0;//not used by BatchLearner
  adaptConfig.threshold_ = 0.001;

  // Set solver during refinement
  solverConfig.eps_ = 0;
  solverConfig.maxIterations_ = 20;
  solverConfig.threshold_ = -1.0;//not used by BatchLearner
  solverConfig.type_ = SGPP::solver::CG;

  // Set parameters for the batchLearner
  batchConfig.filename = "../tests/data/friedman_4d_2000.arff";
  batchConfig.batchsize = 500;
  batchConfig.samples = 500;
  batchConfig.seed = 42;
  batchConfig.wMode = 5;
  batchConfig.wArgument = 1.0;
  batchConfig.refineEvery = 0;
  batchConfig.verbose = true;
  batchConfig.stack = 0;
  batchConfig.testsize = 200;
  batchConfig.lambda = 0.0001f;

  //set up the grid config
  gridConfig.level_ = 4;

  //init the learner
  SGPP::datadriven::BatchLearner learner(batchConfig, gridConfig, solverConfig, adaptConfig);


  int i = 0;

  while (!learner.getIsFinished()) { //while there are still batches to read
    std::cout << "called trainBatch() #" << i << std::endl;
    long int start_l = clock();
    learner.trainBatch();
    long int stop_l = clock();
    cout << "trainBatch() took: " << double(stop_l - start_l) / double(CLOCKS_PER_SEC) * 1000 << endl;

    //If testsize > 0, the learner will test on testsize items following the batch; in the next batch these items will be included and the batch is filled with new items until batchsize is reached. This will be done after every batch.
    //AC = accuracy of the current batch
    //AG = accuracy over all predictions
    cout << "AC: " << i << ", " << learner.getAccCurrent() << "\nAG: " << i << ", " << learner.getAccGlobal() << endl;
    i++;
  }

  //don't forget to close the file (only needed if the file has not been read to completion yet, otherwise the reader will close the file automatically)
  learner.closeFile();

  return EXIT_SUCCESS;
}
