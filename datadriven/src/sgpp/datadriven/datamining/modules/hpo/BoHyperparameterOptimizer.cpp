// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>

#include <vector>
#include <string>
#include <limits>
#include <random>

namespace sgpp {
namespace datadriven {

BoHyperparameterOptimizer::BoHyperparameterOptimizer(SparseGridMiner* miner,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
        : HyperparameterOptimizer(miner, fitterFactory, parser) {
}

double BoHyperparameterOptimizer::run(bool writeToFile) {
  // mute auxiliary optimizers
  base::Printer::getInstance().disableStatusPrinting();

  std::cout << std::endl << "Starting Hyperparameter Optimization using Bayesian Optimization."
      " Results are saved with timestamp." << std::endl << std::endl;

  BOConfig prototype = fitterFactory->getBOConfig();


  // output initialization
  std::ofstream myfile;
  std::stringstream fn;

  if (writeToFile) {
    time_t now = time(nullptr);
    tm tmobj{};
    tm *ltm = localtime_r(&now, &tmobj);
    fn << "Bayesian_" << (ltm->tm_year + 1900) << "_" << (ltm->tm_mon + 1) << "_" << ltm->tm_mday
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
  int bestscnt = 0;
  std::string bestconfigstring;



  // list/vector of configs, start setup
  std::vector<BOConfig> initialConfigs{};
  initialConfigs.reserve(static_cast<size_t>(config.getNRandom()));
  std::mt19937 generator(static_cast<std::mt19937::result_type>(config.getSeed()));

  // random warmup phase
  for (int i = 0; i < config.getNRandom(); ++i) {
    initialConfigs.emplace_back(prototype);
    initialConfigs[i].randomize(generator);
    fitterFactory->setBO(initialConfigs[i]);
    std::string configString = fitterFactory->printConfig();
    // std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    miner->setModel(fitterFactory->buildFitter());
    double result = miner->learn(false);
    // hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    initialConfigs[i].setScore(transformScore(result));
    std::cout << (i + 1) << configString << ", " << result;
    if (writeToFile) {
      myfile.open(fn.str(), std::ios_base::app);
      if (myfile.is_open()) {
        myfile << (i + 1) << configString << ", " << result << std::endl;
      }
      myfile.close();
    }
    if (result < best) {
      best = result;
      bestscnt = i + 1;
      bestconfigstring = configString;
      std::cout << " new best!";
    }
    std::cout << std::endl;
  }

  std::cout << "############# Random Phase finished! #############" << std::endl;

  BayesianOptimization bo(initialConfigs);
  bo.setScales(bo.fitScales(), 0.7);


  // main loop
  for (int q = 0; q < config.getNRuns(); q++) {
    BOConfig nextConfig = bo.main(prototype);
    fitterFactory->setBO(nextConfig);
    std::string configString = fitterFactory->printConfig();
    miner->setModel(fitterFactory->buildFitter());
    double result = miner->learn(false);
    nextConfig.setScore(transformScore(result));
    bo.updateGP(nextConfig, true);
    bo.setScales(bo.fitScales(), 0.1);
    std::cout << (q + config.getNRandom() + 1) << configString << ", " << result;
    if (writeToFile) {
      myfile.open(fn.str(), std::ios_base::app);
      if (myfile.is_open()) {
        myfile << (q + config.getNRandom() + 1) << configString << ", " << result << std::endl;
      }
      myfile.close();
    }
    if (result < best) {
      best = result;
      bestscnt = static_cast<int>(q + config.getNRandom() + 1);
      bestconfigstring = configString;
      std::cout << " new best!";
    }
    std::cout << std::endl;
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

double BoHyperparameterOptimizer::transformScore(double original) {
  return -1 / (1 + original);
}


} /* namespace datadriven */
} /* namespace sgpp */
