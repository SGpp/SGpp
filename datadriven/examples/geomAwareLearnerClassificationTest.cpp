// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_learnerClassificationTest_cpp Learner Classification Test
 * 
 * This represents a small example how to use sparse grids for classification
 * problems. It uses the artificial Ripley dataset.
 */

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/ClassificationLearner.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>
#include <vector>
#include <exception>
#include <limits>
#include <ostream>


std::vector<std::vector<size_t>> getAllInteractions(size_t dimension) {
	size_t geodim = (size_t) std::sqrt(dimension);
	std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();

	for(size_t i = 0; i < geodim-1; i++){
		for(size_t j = 0; j < geodim-1; j++){
			std::vector<size_t> xdir = std::vector<size_t>();
			std::vector<size_t> ydir = std::vector<size_t>();

			xdir.push_back(i*geodim+j);
			ydir.push_back(i*geodim+j);

			xdir.push_back(i*geodim+j+1);
			ydir.push_back((i+1)*geodim+j);

			vec.push_back(xdir);
			vec.push_back(ydir);
		}
	}

	return vec;

}

/**
 * @brief getLearner creates a sparse grid classification learner.
 * @param dimension the number of dimensions.
 * @return a classification for dimension dimension
 */
sgpp::datadriven::ClassificationLearner getLearner(size_t dimension) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  auto adaptivityConfig = sgpp::base::AdpativityConfiguration();
  adaptivityConfig.noPoints_ = 100;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 500;
  solverConfig.eps_ = 1e-8;

  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Diagonal;
  regularizationConfig.lambda_ = 0.00001;
  regularizationConfig.exponentBase_ = 0.25;

  return sgpp::datadriven::ClassificationLearner(gridConfig, adaptivityConfig, solverConfig,
                                                 solverConfig, regularizationConfig, getAllInteractions(dimension));
  //return sgpp::datadriven::ClassificationLearner(gridConfig, adaptivityConfig, solverConfig,
   //                                              solverConfig, regularizationConfig);

}

/**
 * @brief main
   Creates a sparse grid classification learner and prints the training accuracy for the ripley dataset.
 * @return
 */
int main(int argc, char** argv) {
  const auto filenameTrain =
      std::string("../../datasets/house_numbers/SVHN_8x8.extra.arff");


  auto dataTrain = sgpp::datadriven::ARFFTools::readARFF(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();


  std::cout << "Dimensions " << dimensions << "." << std::endl;

  auto learner = getLearner(dimensions);
  learner.train(xTrain, yTrain);
  const auto accuracy = learner.getAccuracy(xTrain, yTrain);
  std::cout << "Best config got a training acc of " << accuracy << "!" << std::endl;

  const auto filenameTest =
      std::string("../../datasets/house_numbers/SVHN_8x8.test.arff");



  auto dataTest = sgpp::datadriven::ARFFTools::readARFF(filenameTest);
  std::cout << "Read file " << filenameTest << "." << std::endl;
  auto xTest = dataTest.getData();
  auto yTest = dataTest.getTargets();
  const auto dimensionsTest = dataTest.getDimension();


  std::cout << "Dimensions " << dimensionsTest << "." << std::endl;

  const auto testAcc = learner.getAccuracy(xTest, yTest);
  std::cout << "Test accuracy of " << testAcc << "!" << std::endl;

  auto predictTest = learner.predict(xTest);

  size_t confusion [10][10] = {0};

  for(size_t i = 0; i < dataTest.getNumberInstances(); i++){
    confusion[(size_t)predictTest.get(i)][(size_t)yTest.get(i)]++;
  }
  std::cout << "Confusion matrix (true\\predicted): " << std::endl;
  for(size_t real = 0; real < 10; real++){
    for(size_t pred = 0; pred < 10; pred++){
      std::cout << confusion[pred][real] << "\t";
    }
    std::cout << std::endl;
  }
}
