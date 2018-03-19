/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HPOScorer.cpp
 *
 *  Created on:	10.12.2017
 *      Author: Eric Koepke
 */
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorer.hpp>

#include <vector>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>

namespace sgpp {
namespace datadriven {

using base::DataVector;

HPOScorer::HPOScorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed,
    double trainPortion, Dataset* testDataset)
    : Scorer{metric, shuffling, seed}, trainPortion{trainPortion}, testDataset{testDataset} {}

// EDIT: disabled because of unique_ptr copying
Scorer* HPOScorer::clone() const {
	// return new HPOScorer{*this};
}

void HPOScorer::resizeTrainData(Dataset& original, Dataset& smaller){
  // perform randomization of indices
  shuffling->setSeed(shuffling->getSeed());
  std::vector<size_t> randomizedIndices(original.getNumberInstances());
  randomizeIndices(original, randomizedIndices);
  Dataset dummyDataset{original.getNumberInstances()-smaller.getNumberInstances(), original.getDimension()};
  splitSet(original, smaller, dummyDataset, randomizedIndices);
}


Dataset* HPOScorer::prepareTestData(Dataset& dataset){
	  // perform randomization of indices
	  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
	  randomizeIndices(dataset, randomizedIndices);
	  // calculate size of testing and training portions
	  size_t trainSize = std::lround(static_cast<double>(dataset.getNumberInstances()) * trainPortion);
	  size_t testSize =  dataset.getNumberInstances() - trainSize;
	  size_t dim = dataset.getDimension();

	  //std::cout<< "Test 1" << std::endl;

	  //Dataset dummyDataset{dataset.getNumberInstances(), dim};

	  // create test and train datasets.
	  //Dataset testDataset{testSize, dim};


	  testDataset = std::make_unique<Dataset>(testSize, dim);
	  Dataset* trainDataset = new Dataset{trainSize, dim};
	  splitSet(dataset, *trainDataset, *testDataset, randomizedIndices);

      //overwrite testdataset with outside test data
      DataSourceBuilder dsbuilder;
      testDataset = std::unique_ptr<Dataset>(dsbuilder.withPath("C:/Users/Eric/Downloads/DE4dTest.csv").assemble()->getNextSamples());

	  return trainDataset;
	  //splitSet(dataset, dummyDataset, trainDataset, randomizedIndices, 2000);
}

void HPOScorer::createTestFile(Dataset& dataset){
  std::vector<size_t> randomizedIndices(dataset.getNumberInstances());
  randomizeIndices(dataset, randomizedIndices);

  // calculate size of testing and training portions
  size_t trainSize = std::lround(static_cast<double>(dataset.getNumberInstances()) * trainPortion);
  size_t testSize =  dataset.getNumberInstances() - trainSize;
  size_t dim = dataset.getDimension();

  testDataset = std::make_unique<Dataset>(testSize, dim);
  Dataset* trainDataset = new Dataset{trainSize, dim};
  splitSet(dataset, *trainDataset, *testDataset, randomizedIndices);

  DataMatrix data(trainDataset->getData());
  DataVector targets(trainDataset->getTargets());
  std::ofstream myfile("C:/Users/Eric/Documents/friedmantrain.txt", std::ios_base::app);
  if(myfile.is_open()) {
    //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;

    for (int i = 0; i < data.getNrows(); i++) {

      for (int k = 0; k < data.getNcols(); k++) {
        myfile << data.get(i,k) << ",";
      }

      myfile << targets[i]  << std::endl;
    }
  }
  myfile.close();
  data = testDataset->getData();
  targets = (testDataset->getTargets());
  myfile = std::ofstream("C:/Users/Eric/Documents/friedmantest.txt", std::ios_base::app);
  if(myfile.is_open()) {
    //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;

    for (int i = 0; i < data.getNrows(); i++) {

      for (int k = 0; k < data.getNcols(); k++) {
        myfile << data.get(i,k) << ",";
      }

      myfile << targets[i]  << std::endl;
    }
  }
  myfile.close();
}

// TODO(lettrich) :recycle
double HPOScorer::calculateScore(ModelFittingBase& model, Dataset& trainDataset,
                                       double* stdDeviation) {

  
  // new code
  // std::cout<< "Test 1" << std::endl;

  // FitterConfiguration* fitterConfig1;
  // fitterConfig1 = fitterFactory->buildConfig(); //*parser
  // ModelFittingBase* model = fitterFactory->buildFitter(fitterConfig1);
  
  // std::cout<< "Test 2" << std::endl;

  
  bool resetVerbose = model.verboseSolver;
  model.verboseSolver = false;
  // int64_t seed = shuffling->getSeed();
  // auto fitterconfig = model->getFitterConfiguration();
  // auto gridConfig = fitterconfig->getGridConfig();
  // auto adaptivityConfig = fitterconfig->getRefinementConfig();
  // auto regularizationConfig = fitterconfig->getRegularizationConfig();
  
  //datadriven::RegularizationConfiguration regularizationConfig;
  //base::AdpativityConfiguration adaptivityConfig;
  //base::RegularGridConfiguration gridConfig;

  // double scores[5][6][5+1][6][9]; //adaptivityConfig.numRefinements_
  // double nscores[4096];
  
  //auto matrix = sgpp::base::DataMatrix(11,0);
  //matrix.appendCol(sgpp::base::DataVector(std::vector<double>({0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1})));
  //auto results = sgpp::base::DataMatrix(11,0);
  //auto temp = sgpp::base::DataVector(std::vector<double>({0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}));

  
  //std::cout<< "Test 2" << std::endl;

/*

  for (int i=1; i <= 4; i++){
    for (int k=1; k <= 5; k++){
      for (int p=1; p <= 5; p++){
        for (int r=1;r <= 8; r++){
          //shuffling->setSeed(seed);
          gridConfig.level_ = i;
          adaptivityConfig.noPoints_ = k;
          adaptivityConfig.threshold_ =  0.0005*p;
          regularizationConfig.lambda_ = pow(10,-r/2.0);
          fitterconfig->setRegularizationConfig( regularizationConfig);
          fitterconfig->setGridConfig( gridConfig);
          fitterconfig->setRefinementConfig( adaptivityConfig);
          model->fit(trainDataset);
          scores[i][k][0][p][r] = test(*model, testDataset);
          std::cout<< "Score: " <<scores[i][k][0][p][r] <<" par: "<<i<<k<<0<<p<<r<< std::endl;
          for (int m=1;m<=5;m++){
            if(model->refine()){
              scores[i][k][m][p][r] = test(*model, testDataset);
              std::cout<< "Score: " <<scores[i][k][m][p][r] <<" par: "<<i<<k<<m<<p<<r<< std::endl;
            }else{
              scores[i][k][m][p][r] = 1;
            }
          }
        }
      }
    }
  }
  */

  /*
  for (int i = 0; i <= 511; i++){
    gridConfig.level_ = (i&3)+1;
    adaptivityConfig.noPoints_ = ((i>>2)&3)+1;
    adaptivityConfig.threshold_ =  0.0005*(((i>>4)&3)+1);
    regularizationConfig.lambda_ = pow(10,-(((i>>6)&7)+1)/2.0);
    fitterconfig->setRegularizationConfig( regularizationConfig);
    fitterconfig->setGridConfig( gridConfig);
    fitterconfig->setRefinementConfig( adaptivityConfig);
    model->fit(trainDataset);
    nscores[i<<3] = test(*model, testDataset);
    std::cout<< "Score: " << nscores[i<<3] <<" par: "<<(i<<3)<< std::endl;
    for (int m=1;m<=7;m++){
      model->refine();
      nscores[(i<<3)+m] = test(*model, testDataset);
      std::cout<< "Score: " <<nscores[(i<<3)+m] <<" par: "<<((i<<3)+m)<< std::endl;
    }
  }
  
  std::cout<< "###################################" << std::endl;
  
  std::ofstream myfile("C:/Users/Eric/Documents/HPOut.csv");
  if(myfile.is_open()){
    myfile<<"refine3,level2,nopoints2,threshold2,lambda3"<<std::endl;
    for (int i = 0; i<=4095; i++){
      for (int k = 0; k <= 11; k++){
        myfile<<(((i>>k)&1)*2-1)<<",";
      }
      myfile<<nscores[i]<<","<<i<<std::endl;
      //if(i>=4092){
      //  break;
      //}
    }
  }
  myfile.close();
  std::cout<< "Write to file finished." << std::endl;
  */


  /*
  std::ofstream myfile("C:/Users/Eric/Documents/HPOut.csv");
  if(myfile.is_open()){
    myfile<<"threshold,noPoints,min,mean,max"<<std::endl;
    for (int p=1;p<=5;p++){
      for (int k=1;k<=5;k++){
        double min = 1;
        double max = 0;
        double sum = 0;
        int cnt = 0;
        for (int m=0;m<=3;m++){
          for (int i=1;i<=4;i++){
            for (int r=1;r<=8;r++){
              if (scores[i][k][m][p][r]<min){
                min = scores[i][k][m][p][r];
              }
              if (scores[i][k][m][p][r]>max){
                max = scores[i][k][m][p][r];
              }
              sum += scores[i][k][m][p][r];
              cnt++;
            }
          }
        }
        myfile<<p<<","<<k<<","<<min<<","<<sum/cnt<<","<<max<<std::endl;
      }
    }
  }
  myfile.close();
  
  double minScore = 1000.0;
   for (int i=1;i<=4;i++){
    for (int k=1;k<=5;k++){
      for (int m=0;m<=3;m++){
        for (int p=1;p<=5;p++){
          for (int r=1;r<=8;r++){
            if (scores[i][k][m][p][r]<=minScore){
              minScore = scores[i][k][m][p][r];
              std::cout<< "Score: " <<minScore <<" par: "<<i<<k<<m<<p<<r<< std::endl;
              gridConfig.level_ = i;
              adaptivityConfig.noPoints_ = k;
              adaptivityConfig.numRefinements_ = m;
              adaptivityConfig.threshold_ =  0.0005*p;
              regularizationConfig.lambda_ = pow(10,-r/2.0);
            }
          }
        }
      }
    }
  }
  
  fitterconfig->setRegularizationConfig( regularizationConfig);
  fitterconfig->setGridConfig( gridConfig);
  fitterconfig->setRefinementConfig( adaptivityConfig);
  */
  
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
  
  
  
  model.fit(trainDataset);
  double score = test(model, *testDataset);
  double best = score+1;
  while(score<best){
	  best = score;
	  model.refine();
	  score = test(model, *testDataset);
	  //std::cout<<"RefinedScore :"<<score<<std::endl;
  }

  model.verboseSolver = resetVerbose;

  if (stdDeviation) {
    *stdDeviation = 0;
  }
  
  // double minScore = 1000.0;

  return best;
}



} /* namespace datadriven */
} /* namespace sgpp */
