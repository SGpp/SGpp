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
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorerFactory.hpp>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource *dataSource,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
    : fitterFactory(fitterFactory) {
  HPOScorerFactory scorerFactory;
  hpoScorer.reset(dynamic_cast<HPOScorer *>(scorerFactory.buildScorer(parser)));
  config.setupDefaults();
  parser.getHPOConfig(config);
  std::unique_ptr<DataSource> ds(dataSource);
  trainData.reset(ds->getNextSamples());
  if (!parser.hasScorerTestset()) {
    trainData.reset(hpoScorer->prepareTestData(*trainData));
  }
  if (config.getNTrainSamples() > 0
      && static_cast<size_t>(config.getNTrainSamples()) <= trainData->getNumberInstances()) {
    Dataset *resize =
        new Dataset(static_cast<size_t>(config.getNTrainSamples()), trainData->getDimension());
    hpoScorer->resizeTrainData(*trainData, *resize);
    trainData.reset(resize);
  }
}

void HyperparameterOptimizer::runHarmonica() {
  Harmonica harmonica{fitterFactory.get()};

  std::cout << std::endl << "Starting Hyperparameter Optimization using Harmonica. Results"
      " are saved with timestamp." << std::endl << std::endl;

  double stdDeviation;   // dummy

  // output initialization
  time_t now = time(nullptr);
  tm tmobj{};  // EDIT: working?
  tm *ltm = localtime_r(&now, &tmobj);
  std::stringstream fn;
  fn << "Harmonica_" << (ltm->tm_year + 1900) << "_" << (ltm->tm_mon + 1) << "_" << ltm->tm_mday
     << "_"
     << ltm->tm_hour << "_" << ltm->tm_min;
  std::ofstream myfile(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
  }
  std::cout << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
  myfile.close();
  double best = std::numeric_limits<double>::infinity();
  int scnt = 1;
  int bestscnt = 0;
  std::string bestconfigstring;

  // loop over stages
  for (size_t q = 0; q < config.getStages().size(); q++) {
    size_t nRuns = (size_t) config.getStages()[q];
    std::vector<std::unique_ptr<ModelFittingBase>> fitters(nRuns);
    DataVector scores(nRuns);
    DataVector transformedScores(nRuns);
    std::vector<std::string> configStrings(nRuns);
    harmonica.prepareConfigs(fitters, static_cast<int>(config.getSeed()), configStrings);

    // run samples (parallelize here)
    for (size_t i = 0; i < nRuns; i++) {
      scores[i] = hpoScorer->calculateScore(*(fitters[i]), *trainData, &stdDeviation);
      std::cout << scnt << configStrings[i] << ", " << scores[i];
      if (scores[i] < best) {
        best = scores[i];
        bestscnt = scnt;
        bestconfigstring = configStrings[i];
        std::cout << " new best!";
      }
      std::cout << std::endl;
      myfile.open(fn.str(), std::ios_base::app);
      if (myfile.is_open()) {
        myfile << scnt << configStrings[i] << ", " << scores[i] << std::endl;
      }
      myfile.close();
      scnt++;
    }

    // constraint introduction
    if (q < config.getStages().size() - 1) {
      harmonica.transformScores(scores, transformedScores);

      harmonica.calculateConstrainedSpace(transformedScores, config.getLambda(),
                                          static_cast<int>(config.getConstraints()[q]));
    }
  }
  myfile.open(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "################## Best Result ###################" << std::endl;
    myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
    myfile << bestscnt << bestconfigstring << ", " << best << std::endl;
    myfile << "##################################################" << std::endl;
  }
  myfile.close();
}

void HyperparameterOptimizer::runBO() {
  // mute auxiliary optimizers
  optimization::Printer::getInstance().disableStatusPrinting();

  std::cout << std::endl << "Starting Hyperparameter Optimization using Bayesian Optimization."
      " Results are saved with timestamp." << std::endl << std::endl;

  BOConfig prototype = fitterFactory->getBOConfig();

  double stdDeviation;   // dummy

  // output initialization
  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::stringstream fn;
  fn << "Bayesian_" << (ltm->tm_year + 1900) << "_" << (ltm->tm_mon + 1) << "_" << ltm->tm_mday
     << "_"
     << ltm->tm_hour << "_" << ltm->tm_min;
  std::ofstream myfile(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
  }
  std::cout << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
  myfile.close();
  double best = std::numeric_limits<double>::infinity();
  int bestscnt = 0;
  std::string bestconfigstring;



  // list/vector of configs, start setup
  std::vector<BOConfig> initialConfigs{};
  initialConfigs.reserve(static_cast<size_t>(config.getNRandom()));
  std::mt19937 generator(static_cast<size_t>(config.getSeed()));

  // random warmup phase
  for (int i = 0; i < config.getNRandom(); ++i) {
    initialConfigs.emplace_back(prototype);// = BOConfig(prototype);
    initialConfigs[i].randomize(generator);
    fitterFactory->setBO(&initialConfigs[i]);
    std::string configString = fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    initialConfigs[i].setScore(result);
    std::cout << (i + 1) << configString << ", " << result;
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << (i + 1) << configString << ", " << result << std::endl;
    }
    myfile.close();
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

  // main loop
  for (int q = 0; q < config.getNRuns(); q++) {
    BOConfig *nextConfig = bo.main(prototype);
    fitterFactory->setBO(nextConfig);
    std::string configString = fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    nextConfig->setScore(result);
    bo.updateGP();
    bo.fitScales();
    std::cout << (q + config.getNRandom() + 1) << configString << ", " << result;
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << (q + config.getNRandom() + 1) << configString << ", " << result << std::endl;
    }
    myfile.close();
    if (result < best) {
      best = result;
      bestscnt = static_cast<int>(q + config.getNRandom() + 1);
      bestconfigstring = configString;
      std::cout << " new best!";
    }
    std::cout << std::endl;
  }
  myfile.open(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "################## Best Result ###################" << std::endl;
    myfile << "SampleNo" << fitterFactory->printHeadline() << ", Score" << std::endl;
    myfile << bestscnt << bestconfigstring << ", " << best << std::endl;
    myfile << "##################################################" << std::endl;
  }
  myfile.close();
}
} /* namespace datadriven */
} /* namespace sgpp */
