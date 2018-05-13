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

HyperparameterOptimizer::HyperparameterOptimizer(DataSource* dataSource, FitterFactory* fitterFactory, DataMiningConfigParser& parser)
    : dataSource(dataSource), fitterFactory(fitterFactory) {
  HPOScorerFactory scorerFactory;
  hpoScorer.reset(static_cast<HPOScorer*>(scorerFactory.buildScorer(parser)));
  config.setupDefaults();
  parser.getHPOConfig(config);
}

void HyperparameterOptimizer::runHarmonica(){
  Harmonica harmonica{fitterFactory.get()}; //EDIT: use correct constructor
  //prepare data
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset dataone{30000, dataset->getDimension()}; //30000
  Dataset datatwo{50000, dataset->getDimension()}; //200
  hpoScorer->resizeTrainData(*dataset, dataone);
  hpoScorer->resizeTrainData(*dataset, datatwo);
  std::cout << "Shuffle test: "<<dataone.getTargets()[0]<<","<<datatwo.getTargets()[0] << std::endl;
  Dataset* trainData = &datatwo;
  //std::unique_ptr<Dataset> trainData = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));
  /*DataSourceBuilder dsbuilder;
  std::unique_ptr<Dataset> readData(dsbuilder.withPath("C:/Users/Eric/Downloads/DR53krandomsamples.csv").assemble()->getNextSamples());
  DataVector readResults(readData->getTargets());
  int seed = (int)std::time(nullptr);//std::random_device()();
  std::cout<<"Seed: "<<seed<<std::endl;
  std::mt19937 sgenerator(seed);
  std::shuffle(readResults.begin(), readResults.end(), sgenerator);
*/

  double stdDeviation;
  double best = 1.0/0; //infinity



for(int q=0;q<config.getStages().size();q++) {
  std::vector<std::unique_ptr<ModelFittingBase>> fitters(config.getStages()[q]);
  DataVector scores(config.getStages()[q]);
  DataVector transformedScores(config.getStages()[q]);
  std::vector<std::string> configStrings(config.getStages()[q]);
  std::vector<int>* configIDs = harmonica.prepareConfigs(fitters, config.getSeed(), configStrings);

  //run samples (parallelize here)

  for (size_t i = 0; i < config.getStages()[q]; i++) {
    //scores[i] = readResults[i];
    scores[i] = hpoScorer->calculateScore(*(fitters[i]), *trainData, &stdDeviation);
    std::cout << "Config: "<<configIDs->at(i)<< ", " << configStrings[i]<<", Score " << i << ":" << scores[i];
    if (scores[i] < best) {
      best = scores[i];
      std::cout << " best!";
    }
    std::cout << std::endl;
   /* std::ofstream myfile("C:/Users/Eric/Documents/DE50kharm.txt", std::ios_base::app);
    if (myfile.is_open()) {
      //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;
      myfile << configIDs->at(i) << ", " << configStrings[i] << ", " << scores[i] << std::endl;
    }
    myfile.close(); */
  }


if(q<config.getStages().size()-1) {


  harmonica.transformScores(scores, transformedScores);

  harmonica.calculateConstrainedSpace(transformedScores, config.getLambda(),
                                      config.getConstraints()[q]); //EDIT: shrink 2 and lambda, normalization
}
  /*std::ofstream myfile("C:/Users/Eric/Documents/dr5harmonica10k.txt", std::ios_base::app);
  if (myfile.is_open()) {
    //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;
    myfile << "#####################################################################" << std::endl;
  }
  myfile.close();*/
}

}



void HyperparameterOptimizer::runBO(){
  optimization::Printer::getInstance().disableStatusPrinting();

  //prepare data
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset dataone{30000, dataset->getDimension()};
  Dataset datatwo{50000, dataset->getDimension()};
  hpoScorer->resizeTrainData(*dataset, dataone);
  hpoScorer->resizeTrainData(*dataset, datatwo);
  std::cout << "Shuffle test: "<<dataone.getTargets()[0]<<","<<datatwo.getTargets()[0] << std::endl;
  Dataset* trainData = &datatwo;

  BOConfig prototype = fitterFactory->getBOConfig();

  double stdDeviation;
  //list/vector of configs, start setup
  std::vector<BOConfig> initialConfigs{};
  initialConfigs.reserve(config.getNRandom());

  std::mt19937 generator(config.getSeed());

  for (int i = 0; i < config.getNRandom(); ++i) {
    initialConfigs.emplace_back(prototype);// = BOConfig(prototype);
    initialConfigs[i].randomize(generator);
    fitterFactory->setBO(&initialConfigs[i]);
    std::string configString = fitterFactory->printConfig();
    std::cout<<configString<<std::endl;
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    initialConfigs[i].setScore(result);
    std::cout << "Result: " << result << std::endl;
   /* std::ofstream myfile("C:/Users/Eric/Documents/DE50kBO.txt", std::ios_base::app);
    if (myfile.is_open()) {
      myfile << configString << ", " <<result << std::endl;
    }
    myfile.close(); */
  }

  std::cout << "Random finished!" << std::endl;


  BayesianOptimization bo(initialConfigs);
  bo.fitScales();


  for(int q=0; q<config.getNRuns(); q++) {
    BOConfig *nextConfig = bo.main(prototype); //EDIT: give prototype earlier
    fitterFactory->setBO(nextConfig);
    std::string configString = fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    std::cout<<configString<<std::endl;
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    nextConfig->setScore(result);
    bo.updateGP();
    std::cout << "Result: " << result << std::endl;
    bo.fitScales();
  /*  std::ofstream myfile("C:/Users/Eric/Documents/DE50kBO.txt", std::ios_base::app);
    if (myfile.is_open()) {
     myfile << configString << ", " <<result << std::endl;
    }
    myfile.close(); */
  }
  // std::cout << "Acquistion: " << min << std::endl;
}


} /* namespace datadriven */
} /* namespace sgpp */
