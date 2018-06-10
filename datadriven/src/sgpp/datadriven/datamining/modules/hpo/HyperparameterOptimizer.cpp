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
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorerFactory.hpp>

#include <iostream>
#include <random>
#include <cmath>

namespace sgpp {
namespace datadriven {

HyperparameterOptimizer::HyperparameterOptimizer(DataSource *dataSource,
                                                 FitterFactory *fitterFactory,
                                                 DataMiningConfigParser &parser)
    : dataSource(dataSource), fitterFactory(fitterFactory) {
  HPOScorerFactory scorerFactory;
  hpoScorer.reset(static_cast<HPOScorer *>(scorerFactory.buildScorer(parser)));
  config.setupDefaults();
  parser.getHPOConfig(config);
}

void HyperparameterOptimizer::runHarmonica() {
  Harmonica harmonica{fitterFactory.get()}; //EDIT: use correct constructor

  //prepare data EDIT: manage data input, data sizes and sources
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset datatwo{50000, dataset->getDimension()}; //200
  hpoScorer->resizeTrainData(*dataset, datatwo);
  Dataset *trainData = &datatwo;


  double stdDeviation; //dummy

  //output initialization
  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::stringstream fn;
  fn << "HPO_output_"<<ltm->tm_year<<"_"<<ltm->tm_mon<<"_"<<ltm->tm_mday<<"_"<<ltm->tm_hour<<"_"
     << ltm->tm_min;
  std::ofstream myfile(fn.str(), std::ios_base::app);
   if (myfile.is_open()) {
     myfile << "SampleNo" <<fitterFactory->printHeadline() << ", Score"<< std::endl;
   }
   myfile.close();
  double best = 1.0 / 0; //infinity
  int scnt = 1;
  int bestscnt = 0;
  std::string bestconfigstring;

  //loop over stages
  for (int q = 0; q < config.getStages().size(); q++) {
    std::vector<std::unique_ptr<ModelFittingBase>> fitters(config.getStages()[q]);
    DataVector scores(config.getStages()[q]);
    DataVector transformedScores(config.getStages()[q]);
    std::vector<std::string> configStrings(config.getStages()[q]);
    std::vector<int>
        *configIDs = harmonica.prepareConfigs(fitters, config.getSeed(), configStrings);

    //run samples (parallelize here)
    for (size_t i = 0; i < config.getStages()[q]; i++) {
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

    //constraint introduction
    if (q < config.getStages().size() - 1) {

      harmonica.transformScores(scores, transformedScores);

      harmonica.calculateConstrainedSpace(transformedScores, config.getLambda(),
                                          config.getConstraints()[q]); //EDIT: normalization
    }
  }
  myfile.open(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "################## Best Result ###################" << std::endl;
    myfile << "SampleNo" <<fitterFactory->printHeadline() << ", Score"<< std::endl;
    myfile << bestscnt << bestconfigstring << ", " << best << std::endl;
    myfile << "##################################################" << std::endl;
  }
  myfile.close();
}

void HyperparameterOptimizer::runBO() {
  //mute auxiliary optimizers
  optimization::Printer::getInstance().disableStatusPrinting();

  //prepare data
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset datatwo{50000, dataset->getDimension()};
  hpoScorer->resizeTrainData(*dataset, datatwo);
  Dataset *trainData = &datatwo;

  BOConfig prototype = fitterFactory->getBOConfig();

  double stdDeviation; //dummy

  //output initialization
  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::stringstream fn;
  fn << "HPO_output_"<<ltm->tm_year<<"_"<<ltm->tm_mon<<"_"<<ltm->tm_mday<<"_"<<ltm->tm_hour<<"_"
     << ltm->tm_min;
  std::ofstream myfile(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "SampleNo" <<fitterFactory->printHeadline() << ", Score"<< std::endl;
  }
  myfile.close();
  double best = 1.0 / 0; //infinity
  int bestscnt = 0;
  std::string bestconfigstring;



  //list/vector of configs, start setup
  std::vector<BOConfig> initialConfigs{};
  initialConfigs.reserve(config.getNRandom());
  std::mt19937 generator(config.getSeed());

  //random warmup phase
  for (int i = 0; i < config.getNRandom(); ++i) {
    initialConfigs.emplace_back(prototype);// = BOConfig(prototype);
    initialConfigs[i].randomize(generator);
    fitterFactory->setBO(&initialConfigs[i]);
    std::string configString = fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    initialConfigs[i].setScore(result);
    std::cout << (i+1) << configString << ", " << result;
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << (i+1) << configString << ", " << result << std::endl;
    }
    myfile.close();
    if (result < best) {
      best = result;
      bestscnt = i+1;
      bestconfigstring = configString;
      std::cout << " new best!";
    }
    std::cout << std::endl;
  }

  std::cout << "############# Random Phase finished! #############" << std::endl;

  BayesianOptimization bo(initialConfigs);

  //main loop
  for (int q = 0; q < config.getNRuns(); q++) {
    BOConfig *nextConfig = bo.main(prototype);
    fitterFactory->setBO(nextConfig);
    std::string configString = fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    nextConfig->setScore(result);
    bo.updateGP();
    bo.fitScales();
    std::cout << (q+config.getNRandom()+1) << configString << ", " << result;
    myfile.open(fn.str(), std::ios_base::app);
    if (myfile.is_open()) {
      myfile << (q+config.getNRandom()+1) << configString << ", " << result << std::endl;
    }
    myfile.close();
    if (result < best) {
      best = result;
      bestscnt = (q+config.getNRandom()+1);
      bestconfigstring = configString;
      std::cout << " new best!";
    }
    std::cout << std::endl;
  }
  myfile.open(fn.str(), std::ios_base::app);
  if (myfile.is_open()) {
    myfile << "################## Best Result ###################" << std::endl;
    myfile << "SampleNo" <<fitterFactory->printHeadline() << ", Score"<< std::endl;
    myfile << bestscnt << bestconfigstring << ", " << best << std::endl;
    myfile << "##################################################" << std::endl;
  }
  myfile.close();
}

} /* namespace datadriven */
} /* namespace sgpp */