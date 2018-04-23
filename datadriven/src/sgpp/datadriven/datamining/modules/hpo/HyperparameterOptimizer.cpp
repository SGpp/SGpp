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
  int seed = 0;

  double stdDeviation;
  double best = 1.0/0; //infinity

  size_t runs[] = {50,50,50}; //EDIT: metaparameter 200,200,100


for(int q=0;q<3;q++) {
  std::vector<std::unique_ptr<ModelFittingBase>> fitters(runs[q]);
  DataVector scores(runs[q]);
  DataVector transformedScores(runs[q]);
  std::vector<std::string> configStrings(runs[q]);
  std::vector<int>* configIDs = harmonica.prepareConfigs(fitters, seed, configStrings);

  //run samples (parallel)

  for (size_t i = 0; i < runs[q]; i++) {
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

  harmonica.transformScores(scores, transformedScores);

  harmonica.calculateConstrainedSpace(transformedScores,0.1,1); //EDIT: shrink 2 and lambda, normalization

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
  initialConfigs.reserve(10);

  std::mt19937 generator(42); //EDIT: meta seed 420 for most tests
  //std::uniform_real_distribution<double> dis(0.0, 1.0);
  //std::cout << "Random test:"<<dis(generator);
  for (int i = 0; i < 10; ++i) {
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


  for(int q=0; q<140; q++) {
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

void HyperparameterOptimizer::runFromFile() { //EDIT: rework, not working
  DataSourceBuilder dsbuilder;
  std::unique_ptr<Dataset> configDataset = std::unique_ptr<Dataset>(
          dsbuilder.withPath("C:/Users/Eric/Documents/fixedconfigs.csv").assemble()->getNextSamples());
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  std::unique_ptr<Dataset> trainData[10];
  trainData[0] = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));
  for (int i = 1; i < 10; i++) {
    trainData[i] = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*trainData[i - 1]));
  }
  //std::unique_ptr<Dataset> trainData = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));


  int nCont;
  std::vector<int> nOptions{};
  fitterFactory->getBOspace(&nCont, nOptions);
  base::DataVector cont(nCont, 0);
  std::vector<int> disc(nOptions.size(), 0);
  double stdDeviation;


  for (int i = 0; i < configDataset->getNumberInstances(); i++) {
    for (int k = 0; k < configDataset->getDimension(); k++) {
      if (k < nCont) {
        cont[k] = configDataset->getData().get(i, k);
      } else {
        disc[k - nCont] = lround(configDataset->getData().get(i, k));
        //std::cout << disc[k] << std::endl; //result
      }
    }
    //EDIT: rework
    //fitterFactory->setBO(cont, disc);
    fitterFactory->printConfig();

    for (int k = 0; k < 10; k++) {
      double result = hpoScorer->calculateScore(*(fitterFactory->buildFitter()), *(trainData[k]), &stdDeviation);
      std::cout << "Result: " << result << std::endl; //result

      //write results in file
      std::ofstream myfile("C:/Users/Eric/Documents/resultsconfig.txt", std::ios_base::app);
      if (myfile.is_open()) {
        //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;

        for (int i = 0; i < nCont; i++) {
          myfile << (cont)[i] << ",";
        }
        for (int i = 0; i < nOptions.size(); i++) {
          myfile << (disc)[i] << ",";
        }
        myfile << trainData[k]->getNumberInstances() << "," << result << std::endl;
      }
      myfile.close();
      //std::cout<< "Write to file finished." << std::endl;
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
