// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>
#include <random>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/RegularizationConfiguration.hpp"
//#include "sgpp/datadriven/application/GaussianKDE.hpp"
//#include "sgpp/datadriven/DatadrivenOpFactory.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;


int main(int argc, char** argv) {
  //std::string filename = "../tests/data/ripleyGarcke.train.arff";
  std::string filename = "../tests/data/banana.arff";

  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& trainData = trainDataset.getData();

  //normalize - banana
  trainData.normalizeDimension(0);
  trainData.normalizeDimension(1);

  // extract train classes
  sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

  //filename = "../tests/data/ripleyGarcke.test.arff";
  filename = "../tests/data/banana.arff";
  // load test samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& testData = testDataset.getData();

  //normalize - banana
  testData.normalizeDimension(0);
  testData.normalizeDimension(1);

  // extract test classes
  sgpp::base::DataVector& testLabels = testDataset.getTargets();

  // configure grid
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = trainDataset.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  //gridConfig.type_ = sgpp::base::GridType::ModLinear;

  // configure adaptive refinement
  std::cout << "# create adaptive refinement config" << std::endl;
  sgpp::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 5;
  adaptConfig.noPoints_ = 4;

  // configure solver
  std::cout << "# create solver config" << std::endl;
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-14;
  solverConfig.threshold_ = 1e-14;

  // configure regularization
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  //regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

  // configure learner
  std::cout << "# create learner config" << std::endl;
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  crossvalidationConfig.enable_ = false;
  crossvalidationConfig.kfold_ = 5;
  crossvalidationConfig.lambda_ = 3.16228e-06;
  crossvalidationConfig.lambdaStart_ = 1e-1;
  crossvalidationConfig.lambdaEnd_ = 1e-10;
  crossvalidationConfig.lambdaSteps_ = 5;
  crossvalidationConfig.logScale_ = true;
  crossvalidationConfig.shuffle_ = true;
  crossvalidationConfig.seed_ = 1234567;
  crossvalidationConfig.silent_ = false;

  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
  learner.initialize(trainData);
  std::cout << "# start training the learner" << std::endl;
  learner.trainOnline(trainLabels);

  std::cout << "# finished training" << std::endl;

  // test learner
  double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
  std::cout << "Acc (train): " << accTrain << std::endl;
  double accTest = learner.getAccuracy(testData, testLabels, 0.0);
  std::cout << "Acc (test): " << accTest << std::endl;

  /*sgpp::base::DataMatrix covSgde(learner.getDim(), learner.getDim());
  learner.cov(covSgde);
  std::cout << covSgde.toString() << std::endl;

  sgpp::datadriven::GaussianKDE kde(samples);
  sgpp::base::DataMatrix covKDE(kde.getDim(), kde.getDim());
  kde.cov(covKDE);
  std::cout << covKDE.toString() << std::endl;

  sgpp::base::DataVector x(learner.getDim());

  for (size_t i = 0; i < x.getSize(); i++) {
    x[i] = 0.5;
  }

  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << learner.getSurpluses()->getSize() << " -> " << learner.getSurpluses()->sum()
            << std::endl;
  std::cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x) << " = pdf_KDE(x)"
            << std::endl;
  std::cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean() << " = mean_KDE(x)"
            << std::endl;
  std::cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance() << " = var_KDE(x)"
            << std::endl;

  sgpp::base::DataMatrix C(gridConfig.dim_, gridConfig.dim_);
  std::cout << "---------------------- Cov_SGDE ------------------------------" << std::endl;
  learner.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "---------------------- Cov KDE--------------------------------" << std::endl;
  kde.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "------------------------------------------------------" << std::endl;
  // inverse Rosenblatt transformation
  auto opInvRos =
      sgpp::op_factory::createOperationInverseRosenblattTransformation(*learner.getGrid().get());
  sgpp::base::DataMatrix points(12, gridConfig.dim_);
  randu(points);

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;

  sgpp::base::DataMatrix pointsCdf(points.getNrows(), points.getNcols());
  opInvRos->doTransformation(learner.getSurpluses().get(), &points, &pointsCdf);

  points.setAll(0.0);
  auto opRos = sgpp::op_factory::createOperationRosenblattTransformation(*learner.getGrid().get());
  opRos->doTransformation(learner.getSurpluses().get(), &pointsCdf, &points);
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << pointsCdf.toString() << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;*/
}
