// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <random>
#include <string>

#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"

void doAllRefinements(SGPP::base::AdpativityConfiguration& adaptConfig, SGPP::base::Grid& grid,
                      SGPP::base::GridGenerator& gridGen, std::mt19937 mt,
                      std::uniform_real_distribution<double>& dist) {
  SGPP::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor myRefineFunc(alphaRefine, adaptConfig.noPoints_,
                                                      adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }
  }
}

int main(int argc, char** argv) {
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "debugging.arff";
  std::string fileName = "friedman_4d.arff";
  //    std::string fileName = "friedman_10d.arff";
  //  std::string fileName = "DR5_train.arff";
  // std::string fileName = "debugging_small.arff";

  uint32_t level = 10;

  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT);

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();
  std::unique_ptr<SGPP::base::Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  SGPP::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;

  SGPP::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;
  std::cout << "number of data points: " << dataset.getNumberInstances() << std::endl;

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());

  for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
    dataSizeVector[i] = static_cast<double>(i + 1);
  }

  std::cout << "creating operation with unrefined grid" << std::endl;
  std::unique_ptr<SGPP::base::OperationMultipleEval> eval =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration);

  doAllRefinements(adaptConfig, *grid, gridGen, mt, dist);

  std::cout << "number of grid points after refinement: " << gridStorage.getSize() << std::endl;
  std::cout << "grid set up" << std::endl;

  SGPP::base::DataVector alphaResult(gridStorage.getSize());

  std::cout << "preparing operation for refined grid" << std::endl;
  eval->prepare();

  std::cout << "calculating result" << std::endl;
  eval->multTranspose(dataSizeVector, alphaResult);

  std::cout << "duration: " << eval->getDuration() << std::endl;

  std::cout << "calculating comparison values..." << std::endl;

  std::unique_ptr<SGPP::base::OperationMultipleEval> evalCompare =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);

  SGPP::base::DataVector alphaResultCompare(gridStorage.getSize());

  evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

  double mse = 0.0;
  for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
    // std::cout << "mine: " << alphaResult[i] << " ref: " <<
    //    alphaResultCompare[i] << std::endl;
    mse += (alphaResult[i] - alphaResultCompare[i]) * (alphaResult[i] - alphaResultCompare[i]);
  }

  mse = mse / static_cast<double>(alphaResult.getSize());
  std::cout << "mse: " << mse << std::endl;
}
