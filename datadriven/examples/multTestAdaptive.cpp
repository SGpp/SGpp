/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <iostream>

#if USE_OCL == 1
#include <random>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

void doAllRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
                      SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen,
                      SGPP::base::DataVector& alpha, std::mt19937 mt,
                      std::uniform_real_distribution<double>& dist) {
  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor* myRefineFunc = new
    SGPP::base::SurplusRefinementFunctor(&alpha,
                                         adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alpha.getSize();
    alpha.resize(grid.getSize());

    for (size_t j = oldSize; j < alpha.getSize(); j++) {
      alpha[j] = dist(mt);
    }

    delete myRefineFunc;
  }
}

int main(int argc, char** argv) {

  /*    SGPP::base::OCLOperationConfiguration parameters;
   parameters.readFromFile("StreamingOCL.cfg");
   std::cout << "internal precision: " << parameters.get("INTERNAL_PRECISION") << std::endl;*/

  //std::string fileName = "debugging_small.arff";
  std::string fileName = "friedman_4d.arff";
  //std::string fileName = "friedman_10d.arff";
  //std::string fileName = "DR5_train.arff";
  //std::string fileName = "debugging.arff";

  std::cout << "Dataset: " << fileName << std::endl;

  SGPP::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", "true");
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", "false");
  parameters.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", "1");
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", "1");
  parameters.addIDAttr("PLATFORM", "first");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", "1");
  parameters.addIDAttr("MAX_DEVICES", "1");
  parameters.addIDAttr("LOCAL_SIZE", "128");
  parameters.addIDAttr("ADAPTIVE_STREAMING_HARD_LIMIT", "10"); //absolute value
  parameters.addIDAttr("ADAPTIVE_STREAMING_DENSITY", "5"); //In percent

  bool bCompare = true;
  bool bBoth = false;

  uint32_t level = 1;

  SGPP::base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 20;
  adaptConfig.numRefinements_ = 2;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::ADAPTIVE,
    SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT);

  SGPP::datadriven::OperationMultipleEvalConfiguration configurationTwo;

  if (argc >= 2) {
    if (strcmp(argv[1], "streaming") == 0) {
      configuration = SGPP::datadriven::OperationMultipleEvalConfiguration(
                        SGPP::datadriven::OperationMultipleEvalType::STREAMING,
                        SGPP::datadriven::OperationMultipleEvalSubType::OCL);
      std::cout << "EvalType::STREAMING" << std::endl;
    } else if ( strcmp(argv[1], "subspace") == 0 ) {
      configuration = SGPP::datadriven::OperationMultipleEvalConfiguration(
                        SGPP::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
                        SGPP::datadriven::OperationMultipleEvalSubType::COMBINED);
      std::cout << "EvalType::SUBSPACECOMBINED" << std::endl;
    } else if ( strcmp(argv[1], "both") == 0 ) {
      configurationTwo = SGPP::datadriven::OperationMultipleEvalConfiguration(
                           SGPP::datadriven::OperationMultipleEvalType::STREAMING,
                           SGPP::datadriven::OperationMultipleEvalSubType::OCL);
      std::cout << "EvalType::ADAPTIVE_AND_STREAMING" << std::endl;
      bBoth = true;
    }
  }

  if ( argc >= 4 ) {
    level = atoi(argv[2]);
    adaptConfig.numRefinements_ = atoi(argv[3]);

    std::cout << "level: " << level << " numRef: " << adaptConfig.numRefinements_ <<
              std::endl;
  }

  if ( argc >= 5) {
    bCompare = static_cast<bool>(atoi(argv[4]));
  }

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();
  SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
  SGPP::base::GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;

  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage->size() << std::endl;
  std::cout << "number of data points: " << dataset.getNumberInstances() <<
            std::endl;

  std::random_device rd;
  std::mt19937 mt(1234);
  std::uniform_real_distribution<double> dist(1, 100);

  SGPP::base::DataVector alpha(gridStorage->size());

  for (size_t i = 0; i < alpha.getSize(); i++) {
    //    alpha[i] = dist(mt);
    alpha[i] = static_cast<double>(i) + 1.0;
  }

  std::cout << "creating operation with unrefined grid" << std::endl;
  SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, trainingData,
        configuration);

  SGPP::base::OperationMultipleEval* evalTwo = nullptr;

  if (bBoth) {
    std::cout << "creating second operation with unrefined grid" << std::endl;
    evalTwo = SGPP::op_factory::createOperationMultipleEval(*grid, trainingData,
              configurationTwo);
  }

  doAllRefinements(adaptConfig, *grid, *gridGen, alpha, mt, dist);

  std::cout << "number of grid points after refinement: " << gridStorage->size()
            << std::endl;
  std::cout << "grid set up" << std::endl;

  SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
  dataSizeVectorResult.setAll(0);

  std::cout << "preparing operation for refined grid" << std::endl;
  eval->prepare();

  std::cout << "calculating result" << std::endl;
  eval->mult(alpha, dataSizeVectorResult);

  if ( bBoth ) {
    std::cout << "preparing second operation for refined grid" << std::endl;
    evalTwo->prepare();
    std::cout << "calculating second result" << std::endl;
    evalTwo->mult(alpha, dataSizeVectorResult);
  }

  std::cout << "durationSub: " << eval->getDuration() << std::endl;

  if ( bBoth ) {
    std::cout << "durationStreaming: " << evalTwo->getDuration() << std::endl;
    std::cout << "diff: " << eval->getDuration() - evalTwo->getDuration() <<
              std::endl;
    std::cout << "speedUp: " << evalTwo->getDuration() / eval->getDuration() <<
              std::endl;
  }

  if ( bCompare) {
    std::cout << "calculating comparison values..." << std::endl;

    SGPP::base::OperationMultipleEval* evalCompare =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);

    SGPP::base::DataVector dataSizeVectorResultCompare(
      dataset.getNumberInstances());
    dataSizeVectorResultCompare.setAll(0.0);

    evalCompare->mult(alpha, dataSizeVectorResultCompare);

    double mse = 0.0;
    double biggestError = 0.0;
    size_t biggestErrorIndex = 0;

    for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
      std::cout << "mine: " << dataSizeVectorResult[i] << " ref: " <<
                dataSizeVectorResultCompare[i] << std::endl;
      double error = (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);

      if (fabs(error) > fabs(biggestError)) {
        biggestError = error;
        biggestErrorIndex = i;
      }

      mse += error * error;
    }

    mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
    std::cout << "mse: " << mse << std::endl;
    std::cout << "biggest error i: " << biggestErrorIndex << " value: " <<
              biggestError << std::endl;
  }

}
#else
int main(int argc, char** argv) {
  std::cout << "need OCL support" << std::endl;
}
#endif
