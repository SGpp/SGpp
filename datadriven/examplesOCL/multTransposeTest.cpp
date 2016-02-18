// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <random>
#include <string>
#include <limits>

#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"

void doAllRefinements(const SGPP::base::AdpativityConfiguration& adaptConfig,
                      SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen, std::mt19937 mt,
                      std::uniform_real_distribution<double>& dist) {
  SGPP::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(
        &alphaRefine, adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }

    delete myRefineFunc;
  }
}

int main(int argc, char** argv) {
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "debugging.arff";
  //  std::string fileName = "friedman_4d.arff";
  std::string fileName = "friedman_10d.arff";
  //  std::string fileName = "DR5_train.arff";
  // std::string fileName = "debugging_small.arff";

  uint32_t level = 7;

  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  SGPP::base::OCLOperationConfiguration parameters("singleDevice.cfg");

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCLMASKMP, parameters);

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  bool modLinear = true;
  SGPP::base::Grid* grid = nullptr;
  if (modLinear) {
    grid = SGPP::base::Grid::createModLinearGrid(dim);
  } else {
    grid = SGPP::base::Grid::createLinearGrid(dim);
  }

  SGPP::base::GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;

  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage->size() << std::endl;
  std::cout << "number of data points: " << dataset.getNumberInstances() << std::endl;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());

  for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
    dataSizeVector[i] = static_cast<double>(i + 1);
  }

  std::cout << "creating operation with unrefined grid" << std::endl;
  SGPP::base::OperationMultipleEval* eval =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration);

  doAllRefinements(adaptConfig, *grid, *gridGen, mt, dist);

  std::cout << "number of grid points after refinement: " << gridStorage->size() << std::endl;
  std::cout << "grid set up" << std::endl;

  SGPP::base::DataVector alphaResult(gridStorage->size());

  std::cout << "preparing operation for refined grid" << std::endl;
  eval->prepare();

  std::cout << "calculating result" << std::endl;

  for (size_t i = 0; i < 10; i++) {
    std::cout << "repeated multTranspose: " << i << std::endl;
    eval->multTranspose(dataSizeVector, alphaResult);
  }

  std::cout << "duration: " << eval->getDuration() << std::endl;

  std::cout << "calculating comparison values..." << std::endl;

  SGPP::base::OperationMultipleEval* evalCompare =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);

  SGPP::base::DataVector alphaResultCompare(gridStorage->size());

  evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

  double mse = 0.0;

  double largestDifferenceMine = 0.0;
  double largestDifferenceReference = 0.0;
  double largestDifference = 0.0;

  for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
    //    std::cout << "mine: " << alphaResult[i] << " ref: " <<
    //    alphaResultCompare[i] << std::endl;
    double difference = std::abs(alphaResult[i] - alphaResultCompare[i]);
    if (difference > largestDifference) {
      largestDifference = difference;
      largestDifferenceMine = alphaResult[i];
      largestDifferenceReference = alphaResultCompare[i];
    }

    std::cout << "difference: " << difference << " mine: " << alphaResult[i]
              << " ref: " << alphaResultCompare[i] << std::endl;

    mse += difference * difference;
  }

  std::cout << "largestDifference: " << largestDifference << " mine: " << largestDifferenceMine
            << " ref: " << largestDifferenceReference << std::endl;

  mse = mse / static_cast<double>(alphaResult.getSize());
  std::cout << "mse: " << mse << std::endl;
}
