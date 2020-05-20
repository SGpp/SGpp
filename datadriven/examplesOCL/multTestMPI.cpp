// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalMPI/OperationMultiEvalMPI.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

#include <mpi.h>

#include <random>
#include <string>

void doAllRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig, sgpp::base::Grid& grid,
                      sgpp::base::GridGenerator& gridGen, sgpp::base::DataVector& alpha) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < adaptivityConfig.numRefinements_; i++) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(alpha, adaptivityConfig.numRefinementPoints_,
                                                      adaptivityConfig.refinementThreshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alpha.getSize();
    alpha.resize(grid.getSize());

    for (size_t j = oldSize; j < alpha.getSize(); j++) {
      alpha[j] = dist(mt);
    }
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*    sgpp::base::OCLOperationConfiguration parameters;
   parameters.readFromFile("StreamingOCL.cfg");
   std::cout << "internal precision: " << parameters.get("INTERNAL_PRECISION")
   << std::endl;*/

  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "debugging.arff";
  std::string fileName = "friedman2_4d_300000.arff";  // NOLINT(build/include_what_you_use)
  //  std::string fileName = "friedman_10d.arff";
  //  std::string fileName = "DR5_train.arff";
  //  std::string fileName = "debugging_small.arff";

  uint32_t level = 10;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  //  sgpp::base::OCLOperationConfiguration parameters("reproduce.cfg");

  //  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
  //      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
  //      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT,
  //      sgpp::datadriven::OperationMultipleEvalMPIType::MASTERSLAVE, parameters);

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT,
      sgpp::datadriven::OperationMultipleEvalMPIType::MASTERSLAVE);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromFile(fileName);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();
  //    std::unique_ptr<sgpp::base::Grid> grid = sgpp::base::Grid::createLinearGrid(dim);

  bool modLinear = false;
  std::unique_ptr<sgpp::base::Grid> grid(nullptr);
  if (modLinear) {
    grid.reset(sgpp::base::Grid::createModLinearGrid(dim));
  } else {
    grid.reset(sgpp::base::Grid::createLinearGrid(dim));
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;
  std::cout << "number of data points: " << dataset.getNumberInstances() << std::endl;

  sgpp::base::DataVector alpha(gridStorage.getSize());

  for (size_t i = 0; i < alpha.getSize(); i++) {
    //    alpha[i] = dist(mt);
    alpha[i] = static_cast<double>(i) + 1.0;
  }

  std::cout << "creating operation with unrefined grid" << std::endl;
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  std::unique_ptr<sgpp::datadriven::OperationMultiEvalMPI> evalMPI(
      dynamic_cast<sgpp::datadriven::OperationMultiEvalMPI*>(eval.release()));

  doAllRefinements(adaptivityConfig, *grid, gridGen, alpha);

  std::cout << "number of grid points after refinement: " << gridStorage.getSize() << std::endl;
  std::cout << "grid set up" << std::endl;

  sgpp::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());

  std::cout << "preparing operation for refined grid" << std::endl;
  evalMPI->prepare();

  std::cout << "calculating result" << std::endl;

  for (size_t i = 0; i < 1; i++) {
    std::cout << "repeated mult: " << i << std::endl;
    if (rank == 0) {
      evalMPI->mult(alpha, dataSizeVectorResult);

      std::cout << "duration: " << evalMPI->getDuration() << std::endl;

      //    sgpp::base::DataVector alpha2(gridStorage.getSize());
      //
      //    eval->multTranspose(dataSizeVectorResult, alpha2);

      std::cout << "calculating comparison values..." << std::endl;

      std::unique_ptr<sgpp::base::OperationMultipleEval> evalCompare(
          sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));

      sgpp::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());

      evalCompare->mult(alpha, dataSizeVectorResultCompare);

      double mse = 0.0;

      double largestDifferenceMine = 0.0;
      double largestDifferenceReference = 0.0;
      double largestDifference = 0.0;

      for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
        double difference = std::abs(dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);
        if (difference > largestDifference) {
          largestDifference = difference;
          largestDifferenceMine = dataSizeVectorResult[i];
          largestDifferenceReference = dataSizeVectorResultCompare[i];
        }

        //    std::cout << "difference: " << difference << " mine: " << dataSizeVectorResult[i]
        //              << " ref: " << dataSizeVectorResultCompare[i] << std::endl;

        mse += difference * difference;
      }

      std::cout << "largestDifference: " << largestDifference << " mine: " << largestDifferenceMine
                << " ref: " << largestDifferenceReference << std::endl;

      mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
      std::cout << "mse: " << mse << std::endl;
    } else {
      evalMPI->multSlave(alpha);
    }
  }

  MPI_Finalize();
}

#else
#include <iostream>
int main(int argc, char** argv) {
  std::cout << "error: build with MPI to enable this example" << std::endl;
  return 0;
}
#endif
