/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * bayesianOptimization.cpp
 *
 *  Created on: June 24, 2018
 *      Author: Eric Koepke
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BOConfig.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>

using sgpp::datadriven::Dataset;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::BOConfig;

BOOST_AUTO_TEST_SUITE(BayesianOptimizationTest)

class ModelFittingTester : public sgpp::datadriven::ModelFittingBase {

  void fit(Dataset &dataset){};

  bool refine(){};

  void update(Dataset &dataset){};

  double evaluate(const DataVector &sample){};

  void evaluate(DataMatrix &samples, DataVector &results){
    results[0] = samples.get(0,0);
  };
};

// EDIT: add FitterConfig and FitterFactory as well to do complete tests
//compare both disctance calculations
//need some concept to test fitScales, sample from a gaussian process maybe?

BOOST_AUTO_TEST_CASE(DistanceCalculation) {
  std::mt19937 generator(123);
  std::vector<int> discOptions = {2,3}; // EDIT: different test case options
  std::vector<int> catOptions = {2,3};
  size_t nCont = 7;
  BOConfig prototype{&discOptions, &catOptions, nCont};
  DataVector scales(prototype.getNPar() + 1, 1);

  double meandistance = 0;
  double mindis = 1000;
  double maxdis = 0;
  BOConfig oldconf(prototype);
  oldconf.randomize(generator);
  for (int j = 0; j < 10000; ++j) {
    BOConfig curconf(prototype);
    curconf.randomize(generator);
    double curdist = oldconf.getScaledDistance(curconf, scales);
    oldconf = curconf;
    meandistance += curdist/10000;
    maxdis = std::fmax(maxdis, curdist);
    mindis = std::fmin(mindis, curdist);
  }
  std::cout << "Mean Distance: " << meandistance << std::endl;
  std::cout << "Min Distance: " << mindis << std::endl;
  std::cout << "Max Distance: " << maxdis << std::endl;
}


BOOST_AUTO_TEST_CASE(addSamples) {
    std::vector<BOConfig> initialConfigs{};
    std::mt19937 generator(123);

    std::vector<int> discOptions = {2,3}; // EDIT: different test case options
    std::vector<int> catOptions = {2,3};
    size_t nCont = 2;
    BOConfig prototype{&discOptions, &catOptions, nCont};

    std::vector<double> scores = {0, 42, 21, 30, 5};
    initialConfigs.reserve(scores.size());


    for (size_t i = 0; i < scores.size(); i++) {
      initialConfigs.emplace_back(prototype);
      initialConfigs[i].randomize(generator);
      initialConfigs[i].setScore(scores[i]);
    }


   //single point not possible because of normalize
   //kernelmatrix both entries necessarry symmetric? yes for tests its necessary
    sgpp::datadriven::BayesianOptimization bo(initialConfigs);

    DataVector scales(initialConfigs.front().getNPar() + 1, 1);

    DataMatrix kernelmatrix(initialConfigs.size(), initialConfigs.size());

  double noise = pow(10, -scales.back() * 10);
  for (size_t i = 0; i < initialConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = bo.kernel(initialConfigs[i].getScaledDistance(initialConfigs[k], scales));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1 + noise);
  }

  DataVector kernelrow(initialConfigs.size());

  DataVector dscores(scores);
  dscores.normalize();
  dscores.sub(DataVector(dscores.size(),
                         dscores.sum() / static_cast<double>(dscores.size())));

    for (size_t i = 0; i < initialConfigs.size(); i++) {
      kernelmatrix.getColumn(i, kernelrow);
      std::cout << "Mean " << i << ": " << bo.mean(kernelrow) << "  |  Original: " << dscores[i] << " | Var: " << bo.var(kernelrow, 1) <<std::endl;
    }

  for (int j = 0; j < 100; ++j) {
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
  }
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
 // std::cout << "kernel value: " << bo.kernel(initialConfigs[0].getScaledDistance(initialConfigs.front(), scales)) << std::endl;
 // std::cout << "weird kernel value: " << initialConfigs[0].getScaledDistance(initialConfigs.front(), scales) << std::endl;
  //  std::cout << "Meany mean: "<< bo.mean(kernelrow) << std::endl;  // EDIT: kself + noise?
}

BOOST_AUTO_TEST_CASE(fitScales) {
  std::vector<BOConfig> initialConfigs{};
  std::mt19937 generator(100);

  std::vector<int> discOptions = {}; // EDIT: different test case options
  std::vector<int> catOptions = {};
  size_t nCont = 2;
  BOConfig prototype{&discOptions, &catOptions, nCont};

  DataVector scores{std::vector<double>({0, 0.5})};

  initialConfigs.reserve(scores.size());

  for (size_t i = 0; i < scores.size(); i++) {
    initialConfigs.emplace_back(prototype);
    initialConfigs[i].randomize(generator);
    initialConfigs[i].setScore(scores[i]);
  }

  //single point not possible because of normalize
  sgpp::datadriven::BayesianOptimization genprocess(initialConfigs);
  sgpp::datadriven::BayesianOptimization fitprocess(initialConfigs);

  DataVector scales{std::vector<double>({1,0.1,0.3})}; //,0.3,0.4,0.01,0.9
  genprocess.setScales(scales, 1);

  for (int j = 0; j < 150; ++j) {
    DataVector kernelrow(initialConfigs.size());
    BOConfig npoint(prototype);
    npoint.randomize(generator);
    for (size_t i = 0; i < initialConfigs.size(); i++) {
      kernelrow[i] = genprocess.kernel(npoint.getScaledDistance(initialConfigs[i], scales));
    }
    double mean = genprocess.mean(kernelrow);
    double var = genprocess.var(kernelrow, 1.001); //EDIT: self variance???
    std::cout << "Mean: " << mean << " | sqVar: " << sqrt(var) << std::endl;
    double delta = scores.max() - scores.min();
    // mean = scores.sum()/scores.size() + mean * delta;
    std::normal_distribution<double> distribution(mean, sqrt(var));
    npoint.setScore(distribution(generator));
    initialConfigs.push_back(npoint);

    scores.push_back(npoint.getScore());
    genprocess.updateGP(npoint);
    fitprocess.updateGP(npoint);
    DataVector fitscales(fitprocess.fitScales());
    double likelihood = fitprocess.likelihood(fitscales);
    //fitscales.componentwise_div(DataVector(fitscales.size(),fitscales.max()));
    std::cout << fitscales.toString() << std::endl;
    fitscales.sub(scales);
    std::cout << "Diff: " << fitscales.l2Norm() << std::endl;
    std::cout << "Likelihood: " << likelihood << " | Original: " << fitprocess.likelihood(scales) << std::endl;
  }
}



BOOST_AUTO_TEST_CASE(validAcquisitionFunction) {
  std::mt19937 generator(123);
  std::uniform_real_distribution<double> ranvar(0.0, 1.0);
  std::uniform_real_distribution<double> ranmean(-1.0, 1.0);
  double oldvar = ranvar(generator);
  double oldmean = ranmean(generator);
  double oldbest = ranmean(generator);
  double oldAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(oldmean, oldvar, oldbest);
  for (int j = 0; j < 1000; ++j) {
    double curvar = ranvar(generator);
    double curmean = ranmean(generator);
    double curbest = ranmean(generator);
    double curAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(curmean, curvar, oldbest);
    if(curmean < oldmean && curvar > oldvar){
      BOOST_CHECK_LT(curAc, oldAc);
      if(oldAc < curAc){
        std::cout << curmean << "," << curvar << "," << curbest << "," << curAc << std::endl;
        std::cout << oldmean << "," << oldvar << "," << oldbest << "," << oldAc << std::endl ;
        std::cout << std::endl;
      }
    }
    oldvar = curvar;
    oldmean = curmean;
    oldbest = curbest;
    oldAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(curmean, curvar, curbest);
  }

  std::cout << sgpp::datadriven::BayesianOptimization::acquisitionEI(-0.5, 0, 0.9) << std::endl;
  std::cout << sgpp::datadriven::BayesianOptimization::acquisitionEI(-0.50001, 0.8, 0.9) << std::endl;
}



BOOST_AUTO_TEST_SUITE_END()
