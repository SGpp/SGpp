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
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/OperationMultipleEvalMatrix.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>


#include <iostream>
#include <random>
#include <cmath>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BayesianOptimization.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include "Harmonica.hpp"
#include "HPOScorerFactory.hpp"

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
  Dataset dataone{10000, dataset->getDimension()}; //30000
  Dataset datatwo{2000, dataset->getDimension()}; //200
  hpoScorer->resizeTrainData(*dataset, dataone);
  hpoScorer->resizeTrainData(*dataset, datatwo);
  std::cout << "Shuffle test: "<<dataone.getTargets()[0]<<","<<datatwo.getTargets()[0] << std::endl;
  Dataset* trainData = &dataone;
  //std::unique_ptr<Dataset> trainData = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));

  double stdDeviation;
  double best = 1.0/0; //infinity

  int n = 200; //EDIT: metaparameter
  std::vector<std::unique_ptr<ModelFittingBase>> fitters(n);
  DataVector scores(n);
  DataVector transformedScores(n);

for(int q=0;q<4;q++) {
  harmonica.prepareConfigs(fitters);

  //run samples (parallel)

  for (int i = 0; i < n; i++) {
    scores[i] = hpoScorer->calculateScore(*(fitters[i]), *trainData, &stdDeviation);
    std::cout << "Score " << i << ":" << scores[i];
    if (scores[i] < best) {
      best = scores[i];
      std::cout << " best!";
    }
    std::cout << std::endl;
  }

  harmonica.transformScores(scores, transformedScores);

  harmonica.calculateConstrainedSpace(transformedScores,0.001,2);
}

}



void HyperparameterOptimizer::runFromFile(){ //EDIT: rework, not working
  DataSourceBuilder dsbuilder;
  std::unique_ptr<Dataset> configDataset = std::unique_ptr<Dataset>(dsbuilder.withPath("C:/Users/Eric/Documents/fixedconfigs.csv").assemble()->getNextSamples());
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  std::unique_ptr<Dataset> trainData[10];
  trainData[0] = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));
  for(int i=1;i<10;i++){
    trainData[i] = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*trainData[i-1]));
  }
  //std::unique_ptr<Dataset> trainData = std::unique_ptr<Dataset>(hpoScorer->prepareTestData(*dataset));


  int nCont;
  std::vector<int> nOptions{};
  fitterFactory->getBOspace(&nCont, nOptions);
  base::DataVector cont(nCont, 0);
  std::vector<int> disc(nOptions.size(), 0);
  double stdDeviation;


  for(int i=0; i<configDataset->getNumberInstances();i++){
    for(int k=0; k<configDataset->getDimension(); k++){
      if(k<nCont) {
        cont[k] = configDataset->getData().get(i, k);
      }else {
        disc[k-nCont] = lround(configDataset->getData().get(i, k));
        //std::cout << disc[k] << std::endl; //result
      }
    }
    //EDIT: rework
    //fitterFactory->setBO(cont, disc);
    fitterFactory->printConfig();

    for(int k=0; k<10; k++){
      double result = hpoScorer->calculateScore(*(fitterFactory->buildFitter()), *(trainData[k]), &stdDeviation);
      std::cout << "Result: " << result << std::endl; //result

      //write results in file
      std::ofstream myfile("C:/Users/Eric/Documents/resultsconfig.txt", std::ios_base::app);
      if(myfile.is_open()) {
        //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;

        for (int i = 0; i < nCont; i++) {
          myfile << (cont)[i] <<",";
        }
        for (int i = 0; i < nOptions.size(); i++) {
          myfile << (disc)[i]<<",";
        }
        myfile << trainData[k]->getNumberInstances() << "," << result << std::endl;
      }
      myfile.close();
      //std::cout<< "Write to file finished." << std::endl;
    }
  }
}

void HyperparameterOptimizer::runBO(){
  optimization::Printer::getInstance().disableStatusPrinting();

  //prepare data
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset dataone{30000, dataset->getDimension()};
  Dataset datatwo{200, dataset->getDimension()};
  hpoScorer->resizeTrainData(*dataset, dataone);
  hpoScorer->resizeTrainData(*dataset, datatwo);
  std::cout << "Shuffle test: "<<dataone.getTargets()[0]<<","<<datatwo.getTargets()[0] << std::endl;
  Dataset* trainData = &dataone;

  BOConfig prototype = fitterFactory->getBOConfig();

  double stdDeviation;
  //list/vector of configs, start setup
  std::vector<BOConfig> initialConfigs{};
  initialConfigs.reserve(10);

  std::mt19937 generator(42); //EDIT: meta seed
  //std::uniform_real_distribution<double> dis(0.0, 1.0);
  //std::cout << "Random test:"<<dis(generator);
  for (int i = 0; i < 10; ++i) {
    initialConfigs.emplace_back(prototype);// = BOConfig(prototype);
    initialConfigs[i].randomize(generator);
    fitterFactory->setBO(&initialConfigs[i]);
    fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    initialConfigs[i].setScore(result);
    std::cout << "Result: " << result << std::endl;
  }

  std::cout << "Random finished!" << std::endl;


  BayesianOptimization bo(initialConfigs);

  for(int q=0; q<500; q++) {
    BOConfig *nextConfig = bo.main(prototype); //EDIT: give prototype earlier
    fitterFactory->setBO(nextConfig);
    fitterFactory->printConfig();
    std::unique_ptr<ModelFittingBase> fitter(fitterFactory->buildFitter());
    double result = hpoScorer->calculateScore(*fitter, *trainData, &stdDeviation);
    nextConfig->setScore(result);
    bo.updateGP();
    std::cout << "Result: " << result << std::endl;
  }
  // std::cout << "Acquistion: " << min << std::endl;
}

void HyperparameterOptimizer::testBO() {
  double kernelwidth = 1.5;

  optimization::Printer::getInstance().disableStatusPrinting();

  //prepare data
  std::unique_ptr<Dataset> dataset(dataSource->getNextSamples());
  Dataset dataone{30000, dataset->getDimension()};
  Dataset datatwo{200, dataset->getDimension()};
  hpoScorer->resizeTrainData(*dataset, dataone);
  hpoScorer->resizeTrainData(*dataset, datatwo);
  std::cout << "Shuffle test: "<<dataone.getTargets()[0]<<","<<datatwo.getTargets()[0] << std::endl;
  Dataset* trainData = &datatwo;

  std::vector<base::DataVector> contPoints{};
  std::vector<std::vector<int>> discPoints{};
  base::DataVector results{};
  double bestsofar = 1000;

  int nCont;
  std::vector<int> nOptions{};
  fitterFactory->getBOspace(&nCont, nOptions);
  int total = 1;
  for(int max: nOptions){
    total = total*max;
  }

  std::cout << "Total:" << total << std::endl;



  base::DataVector cont(nCont, 0);
  std::vector<int> discrete(nOptions.size(), 0);
  //fitterFactory->setBO(cont, discrete);
  double stdDeviation;
  bestsofar = -1/(1+hpoScorer->calculateScore(*(fitterFactory->buildFitter()), *trainData, &stdDeviation)); //-300
  std::vector<BOConfig> initialConfigs{1};
  BayesianOptimization BO(initialConfigs);

  results.push_back(bestsofar);
  contPoints.push_back(cont);
  discPoints.push_back(discrete);

  for(int q=0; q<1000; q++) {


    std::vector<int> disc(nOptions.size(), 0);

    double min = 1.0/0; // inf
    std::unique_ptr<base::DataVector> mincon;
    std::unique_ptr<std::vector<int>> mindisc;
    for (int i = 0; i < total; i++) {
      int k = i;
      int p = 0;
      for (int max: nOptions) {
        disc[p] = k % max; //TEST: does this work?
        k = k / max;
        p++;
      }
      std::vector<double> squaresum(discPoints.size(), 0);
      for (int i = 0; i < discPoints.size(); i++) {
        for (int k = 0; k < nOptions.size(); k++) {
          squaresum[i] += std::pow(static_cast<double>(disc[k] - discPoints[i][k]) / (nOptions[k] - 1), 2);
        }
      }
      std::function<double(const base::DataVector &)> func =
              [&squaresum, &contPoints, &BO, &bestsofar, &kernelwidth](const base::DataVector &inp) {
                  base::DataVector kernels(squaresum.size());
                  for (int i = 0; i < squaresum.size(); i++) {
                    base::DataVector tmp(inp);
                    tmp.sub(contPoints[i]);
                    kernels[i] = exp((0-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2 * kernelwidth); // divided by 2
                    if(kernels[i]==1){
                      return 1.0/0;
                    }
                    // std::cout << "Kernel value: "<<exp((-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2) <<std::endl;
                  }
                  return BO.acquisitionEI(kernels, 1, bestsofar);  //EDIT: ascend or descent?
                  //return BO.var(kernels, 1);
              };
      // std::cout << "Test Point " << i <<", min: "<<min<<std::endl;
      /* bool pronty = true;
      for (int k = 0; k < disc.size(); k++) {
        if(disc[k]!=discPoints[0][k]){
          pronty = false;
        }
      }
      if(pronty){
        std::cout << "Variance of point 0:" << func(contPoints[0]) << std::endl;
        base::DataVector kernels(squaresum.size());
        for (int i = 0; i < squaresum.size(); i++) {
          base::DataVector tmp(contPoints[0]);
          tmp.sub(contPoints[i]);
          kernels[i] = exp((0 - squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2);
          std::cout << "Kernel value: "<<exp((-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2) <<std::endl;
        }
        } */
      optimization::WrapperScalarFunction wrapper(nCont, func);
      optimization::optimizer::MultiStart optimizer(wrapper);
      optimizer.optimize();
      if (optimizer.getOptimalValue() < min) {
        min = optimizer.getOptimalValue();
        mincon = std::make_unique<base::DataVector>(optimizer.getOptimalPoint());
        mindisc = std::make_unique<std::vector<int>>(disc);
      }
      // std::cout<<optimizer.getOptimalPoint()[0]<<","<<optimizer.getOptimalPoint()[1]<<":"<<optimizer.getOptimalValue()<<std::endl;
    }
    //fitterFactory->setBO(*mincon, *mindisc);
    double result = -1/(1+hpoScorer->calculateScore(*(fitterFactory->buildFitter()), *trainData, &stdDeviation)); //-300
    fitterFactory->printConfig();
    std::cout << "Acquistion: " << min << std::endl;
    std::cout << "Result: " << -1/result-1 << std::endl; //result

    //write results in file
    /*
    std::ofstream myfile("C:/Users/Eric/Documents/BOresultsDE.txt", std::ios_base::app);
    if(myfile.is_open()) {
      //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;

      for (int i = 0; i < nCont; i++) {
        myfile << (*mincon)[i] <<",";
      }
      for (int i = 0; i < nOptions.size(); i++) {
        myfile << (*mindisc)[i]<<",";
      }
      myfile << -1/result-1 << std::endl;
    }
    myfile.close();
    std::cout<< "Write to file finished." << std::endl;
    */
    if(result<bestsofar){
      bestsofar = result;
    }

    results.push_back(result);
    contPoints.push_back(*mincon);
    discPoints.push_back(*mindisc);

    std::vector<double> squaresum(discPoints.size(), 0);
    for (int i = 0; i < discPoints.size(); i++) {
      for (int k = 0; k < nOptions.size(); k++) {
        squaresum[i] += std::pow(static_cast<double>((*mindisc)[k] - discPoints[i][k]) / (nOptions[k] - 1), 2);
      }
    }
    base::DataVector kernels(squaresum.size());
    for (int i = 0; i < squaresum.size(); i++) {
      base::DataVector tmp(*mincon);
      tmp.sub(contPoints[i]);
      kernels[i] = exp((-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2 * kernelwidth); //devided by 2
    }
    kernels[kernels.size()-1]=1.0001; //EDIT: noise added
    //BO.updateGP(kernels, results);  //EDIT: ascend or descent?

  }
  std::cout << "Best: " << -1/bestsofar-1 << std::endl; //+300


}

} /* namespace datadriven */
} /* namespace sgpp */
