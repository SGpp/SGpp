/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFactory.cpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp>


//#include <math.h>

#include <string>

namespace sgpp {
namespace datadriven {

SparseGridMiner* LeastSquaresRegressionMinerFactory::buildMiner(const std::string& path) const {
  DataMiningConfigParser parser(path);

  return new SparseGridMiner(createDataSource(parser), createFitter(parser), createScorer(parser));
}

void LeastSquaresRegressionMinerFactory::optimizeHyperparameters(const std::string& path){
  DataMiningConfigParser parser(path);
  
  double scores [5][5][6][6][9];
  //double value1 [5];
  //double value2 [5];
  //double value3 [5];
  auto matrix = sgpp::base::DataMatrix(11,0);
  matrix.appendCol(sgpp::base::DataVector(std::vector<double>({0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1})));
  auto results = sgpp::base::DataMatrix(11,0);
  auto temp = sgpp::base::DataVector(std::vector<double>({0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}));

  
  TestingConfiguration testconfig;
  parser.getScorerTestingConfig(testconfig, testconfig);
 // std::unique_ptr<ScorerFactory> scfactory;
 // scfactory = std::make_unique<SplittingScorerFactory>();
  auto metric = new MSE();
  auto shuffling = new RandomShufflingFunctor();
  Scorer* sc = new SplittingScorer(metric, shuffling, testconfig.randomSeed, testconfig.testingPortion);
  // auto ds = createDataSource(parser);
  // std::unique_ptr<Dataset> dataset(ds->getNextSamples());
  auto dataset = sgpp::datadriven::CSVTools::readCSV("C:/Users/Eric/Documents/Test3.csv", true);
  double stdDeviation;
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  ModelFittingLeastSquares* model;


  for (int i=1; i <= 4; i++){
    for (int k=0; k <= 4; k++){
      for (int m= 1;m <= 5; m++){
        for (int p=1; p <= 5; p++){
          for (int r=1;r <= 8; r++){
            shuffling->setSeed(testconfig.randomSeed);
            config.setHyperParameters(i, k, m, 0.0005*p, pow(10,-r/2.0));
            model = (new ModelFittingLeastSquares(config));
            scores[i][k][m][p][r] = sc->calculateScore(*model, dataset, &stdDeviation);
          }
        }
      }
    }
  }
double minScore = 1.0;
   for(int i=1;i<=4;i++){
    for(int k=0;k<=4;k++){
      for(int m=1;m<=5;m++){
        for(int p=1;p<=5;p++){
          for(int r=1;r<=8;r++){
            if(scores[i][k][m][p][r]<minScore){
              minScore = scores[i][k][m][p][r];
              std::cout<< "Score: " <<minScore <<" par: "<<i<<k<<m<<p<<r<< std::endl;
            }
          }
        }
      }
    }
  }


 /* for(int i=0;i<5;i++){
    config.setRefinements(i);
    model = (new ModelFittingLeastSquares(config));
    scores[i] = sc->calculateScore(*model, dataset, &stdDeviation);
    //value1[i] = model->evaluate(sgpp::base::DataVector(std::vector<double>({0})));
    //value2[i] = model->evaluate(sgpp::base::DataVector(std::vector<double>({0.86})));
    model->evaluate(matrix, temp);
    results.appendCol(temp);
  }
  for(int i=0;i<5;i++){
    std::cout<< "###############" << std::endl
    << "Score: " <<scores[i] << std::endl;
    //<< "Point 0: "<<value1[i] << std::endl
    //<< "Point 0: "<<value1[i] << std::endl
    //<< "Point 0.86: "<<value2[i] << std::endl;
    for(int k=0;k<11;k++){
    std::cout<< "Point"<<k<<": " <<results.get(k,i) << std::endl;
    //<< "Point 0: "<<value1[i] << std::endl
    //<< "Point 0: "<<value1[i] << std::endl
    //<< "Point 0.86: "<<value2[i] << std::endl;
    }
  }
  */
  
  /*
  auto ds = createDataSource(parser);
  auto sc = createScorer(parser);
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  
  config.setRefinements(1);
  ModelFittingLeastSquares fitter(config);
  SparseGridMiner miner(ds, &fitter, sc);
  miner.learn();
  std::cout << "ERSTERFERTIG" << std::endl;
  auto ds2 = createDataSource(parser);
  auto sc2 = createScorer(parser);
  FitterConfigurationLeastSquares config2{};
  config2.readParams(parser);
  config2.setRefinements(4);
  auto fitter2 = ModelFittingLeastSquares(config2);
  miner = SparseGridMiner(ds2, &fitter2, sc2);
  miner.learn();
  */
}

DataSource* LeastSquaresRegressionMinerFactory::createDataSource(
    const DataMiningConfigParser& parser) const {
  DataSourceConfig config;

  bool hasSource = parser.getDataSourceConfig(config, config);

  if (hasSource && config.filePath.compare("") != 0) {
    DataSourceBuilder builder;
    return builder.fromConfig(config);
  } else {
    throw base::data_exception("No file name provided for datasource.");
  }
}

ModelFittingBase* LeastSquaresRegressionMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  return new ModelFittingLeastSquares(config);
}

Scorer* LeastSquaresRegressionMinerFactory::createScorer(
    const DataMiningConfigParser& parser) const {
  std::unique_ptr<ScorerFactory> factory;

  if (parser.hasScorerConfigCrossValidation()) {
    factory = std::make_unique<CrossValidationScorerFactory>();
  } else {
    factory = std::make_unique<SplittingScorerFactory>();
  }
  return factory->buildScorer(parser);
}

} /* namespace datadriven */
} /* namespace sgpp */
