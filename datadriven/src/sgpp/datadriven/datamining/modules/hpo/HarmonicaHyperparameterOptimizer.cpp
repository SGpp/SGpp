// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>

#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>

#include <vector>
#include <string>
#include <limits>

namespace sgpp {
namespace datadriven {

HarmonicaHyperparameterOptimizer::HarmonicaHyperparameterOptimizer(SparseGridMiner* miner,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
        : HyperparameterOptimizer(miner, fitterFactory, parser) {
}


double HarmonicaHyperparameterOptimizer::run(bool writeToFile) {
  Harmonica harmonica{fitterFactory.get()};

  std::cout << std::endl << "Starting Hyperparameter Optimization using Harmonica. Results"
          " are saved with timestamp." << std::endl << std::endl;


  // output initialization
  std::ofstream myfile;
  std::stringstream fn;

  if (writeToFile) {
    time_t now = time(nullptr);
    tm tmobj{};
    tm *ltm = localtime_r(&now, &tmobj);
    fn << "Harmonica_" << (ltm->tm_year + 1900) << "_" << (ltm->tm_mon + 1) << "_" << ltm->tm_mday
       << "_"
       << ltm->tm_hour << "_" << ltm->tm_min;
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
    }
    std::cout << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
    myfile.close();
  }
  double best = std::numeric_limits<double>::infinity();
  int scnt = 1;
  int bestscnt = 0;
  std::string bestconfigstring;

  // loop over stages
  for (size_t q = 0; q < config.getStages().size(); q++) {
    size_t nRuns = static_cast<size_t>(config.getStages()[q]);
    std::vector<ModelFittingBase*> fitters(nRuns);
    DataVector scores(nRuns);
    DataVector transformedScores(nRuns);
    std::vector<std::string> configStrings(nRuns);
    harmonica.prepareConfigs(fitters, static_cast<int>(config.getSeed()), configStrings);

    // run samples (parallelize here)
    for (size_t i = 0; i < nRuns; i++) {
      miner->setModel(fitters[i]);
      scores[i] = miner->learn(false);
      // hpoScorer->calculateScore(*(fitters[i]), *trainData, &stdDeviation);
      std::cout << scnt << configStrings[i] << ", " << scores[i];
      if (scores[i] < best) {
        best = scores[i];
        bestscnt = scnt;
        bestconfigstring = configStrings[i];
        std::cout << " new best!";
      }
      std::cout << std::endl;
      if (writeToFile) {
        myfile.open(fn.str(), std::ios_base::app);
        if (myfile.is_open()) {
          myfile << scnt << configStrings[i] << ", " << scores[i] << std::endl;
        } else {
          std::cout << "Output File '" << fn.str() << "' can't be written to." << std::endl;
        }
        myfile.close();
      }
      scnt++;
    }

    // constraint introduction
    if (q < config.getStages().size() - 1) {
      harmonica.transformScores(scores, transformedScores);

      harmonica.calculateConstrainedSpace(transformedScores, config.getLambda(),
                                          static_cast<int>(config.getConstraints()[q]));
    }
  }
  if (writeToFile) {
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << "################## Best Result ###################" << std::endl;
      myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
      myfile << bestscnt << bestconfigstring << ", " << best << std::endl;
      myfile << "##################################################" << std::endl;
    }
    myfile.close();
  }
  return best;
}
} /* namespace datadriven */
} /* namespace sgpp */
