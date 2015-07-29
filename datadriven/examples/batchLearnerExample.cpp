// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/BatchLearner.hpp>
#include <sgpp/datadriven/application/BatchConfiguration.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

/**
 * This programm demonstrates the usage of the BatchLearner class. After the parameters are set, the method trainBatch() is called until the end of the file has been reached.
 */

using namespace sg::base;
using namespace sg::datadriven;
using namespace std;


int main (int argc, char** args) {
  cout << "parameters: bs(batch size), ts (test size), input, mode(weighting), arg(for weighting), level, pts(to refine every refinement), ref(refine every xth batch, 0=never)" << endl;

  std::map<string, string> argsMap;

  for (int i = 1; i < argc; i += 2) {
    argsMap[args[i]] = args[i + 1];
  }


  //set variables
  sg::base::BatchConfiguration batchConfig;
  sg::solver::SLESolverConfiguration solverConfig;
  sg::base::AdpativityConfiguration adaptConfig;
  sg::base::RegularGridConfiguration gridConfig;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;//not used by BatchLearner
  adaptConfig.noPoints_ = atoi(argsMap["pts"].c_str());
  adaptConfig.numRefinements_ = 1;//not used by BatchLearner
  adaptConfig.percent_ = 10.0;//not used by BatchLearner
  adaptConfig.threshold_ = 0.001;

  // Set solver during refinement
  solverConfig.eps_ = 0;
  solverConfig.maxIterations_ = 20;
  solverConfig.threshold_ = -1.0;//not used by BatchLearner
  solverConfig.type_ = sg::solver::CG;

  // Set parameters for the batchLearner
  batchConfig.filename = argsMap["input"].c_str();
  batchConfig.batchsize = atoi(argsMap["bs"].c_str());
  batchConfig.samples = 500;
  batchConfig.seed = 42;
  batchConfig.wMode = atoi(argsMap["mode"].c_str());;
  batchConfig.wArgument = atof(argsMap["arg"].c_str());
  batchConfig.refineEvery = atoi(argsMap["ref"].c_str());
  batchConfig.verbose = true;
  batchConfig.stack = 0;
  batchConfig.testsize = atoi(argsMap["ts"].c_str());
  batchConfig.lambda = 0.0001f;

  //set up the grid config
  gridConfig.level_ = atoi(argsMap["level"].c_str());

  //init the learner
  sg::datadriven::BatchLearner learner(batchConfig, gridConfig, solverConfig, adaptConfig);


  int i = 0;

  while (!learner.getIsFinished()) { //while there are still batches to read
    std::cout << "called trainBatch() #" << i << std::endl;
    long int start_l = clock();
    learner.trainBatch();
    long int stop_l = clock();
    cout << "trainBatch() took: " << (stop_l - start_l) / double(CLOCKS_PER_SEC) * 1000 << endl;

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
