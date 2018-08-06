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

BOOST_AUTO_TEST_CASE(addSamples) {
    std::vector<BOConfig> initialConfigs{};
    initialConfigs.reserve(1);
    std::mt19937 generator(123);

    std::vector<int> discOptions = {2,3}; // EDIT: different test case options
    std::vector<int> catOptions = {2,3};
    size_t nCont = 2;
    BOConfig prototype{&discOptions, &catOptions, nCont};
    initialConfigs.emplace_back(prototype);
    initialConfigs[0].randomize(generator);
    initialConfigs[0].setScore(42);
    initialConfigs.emplace_back(prototype);
    initialConfigs[1].randomize(generator);
    initialConfigs[1].setScore(0);
    initialConfigs.emplace_back(prototype);
    initialConfigs[2].randomize(generator);
    initialConfigs[2].setScore(21);

   //single point not possible because of normalize
   //kernelmatrix both entries necessarry symmetric?
    sgpp::datadriven::BayesianOptimization bo(initialConfigs);

    DataVector scales(initialConfigs.front().getNPar() + 1, 1);


  DataVector kernelrow(initialConfigs.size());
    for (size_t i = 0; i < initialConfigs.size(); i++) {
      kernelrow[i] = bo.kernel(initialConfigs[0].getScaledDistance(initialConfigs.front(), scales));
    }
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  std::cout << "kernel value: " << bo.kernel(initialConfigs[0].getScaledDistance(initialConfigs.front(), scales)) << std::endl;
  std::cout << "weird kernel value: " << initialConfigs[0].getScaledDistance(initialConfigs.front(), scales) << std::endl;
    std::cout << "Meany mean: "<< bo.mean(kernelrow) << std::endl;  // EDIT: kself + noise?
}

BOOST_AUTO_TEST_SUITE_END()
