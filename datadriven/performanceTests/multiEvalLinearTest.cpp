// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include "testsCommon.hpp"

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

#include <zlib.h>
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#define OUT_FILENAME "results.csv"
// #define REFINEMENT_POINTS 100

// std::vector<std::string> fileNames = { "datadriven/datasets/friedman/friedman2_4d_10000.arff.gz",
//        "datadriven/datasets/friedman/friedman1_10d_2000.arff.gz",
//        "datadriven/datasets/DR5/DR5_train.arff.gz" };
//
// std::vector<std::string> datasetNames = { "Friedman 4d", "Friedman 10d", "DR5" };
//
// std::vector<size_t> levels = { 9, 5, 7 };
// std::vector<size_t> refinementSteps = { 70, 70, 70 };
//
// std::vector<size_t> levelsModLinear = { 9, 5, 7 };
// std::vector<size_t> refinementStepsModLinear = { 70, 70, 70 };

std::vector<std::string> fileNames = {"datadriven/datasets/DR5/DR5_train.arff.gz"};

std::vector<std::string> datasetNames = {"DR5"};

std::vector<size_t> levels = {7};

std::vector<size_t> levelsModLinear = {7};

std::map<std::string, sgpp::datadriven::MetaLearner*> preparedGrids;

std::map<std::string, sgpp::datadriven::MetaLearner*> preparedGridsModLinear;

struct HPCSE2015Fixture {
  HPCSE2015Fixture() {
    BOOST_TEST_MESSAGE("setup fixture");
    outFile.open(OUT_FILENAME);
    outFile << "Dataset, Basis, Kernel, Grid size, Duration (s)" << std::endl;
  }
  ~HPCSE2015Fixture() {
    BOOST_TEST_MESSAGE("teardown fixture");
    outFile.close();
  }
  std::ofstream outFile;
} logger;

static size_t refinedGridSize = 0;

void getRuntime(sgpp::base::GridType gridType, const std::string& kernel, std::string& fileName,
                std::string& datasetName, size_t level,
                sgpp::base::AdaptivityConfiguration adaptivityConfig,
                sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  } else {
    throw;
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage.getDimension());

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);
  BOOST_TEST_MESSAGE("number of grid points: " << gridStorage.getSize());
  BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());

  doDirectedRefinements(adaptivityConfig, *grid, gridGen);
  //    doRandomRefinements(adaptivityConfig, *grid, gridGen);

  BOOST_TEST_MESSAGE("size of refined grid: " << gridStorage.getSize());
  refinedGridSize = gridStorage.getSize();

  sgpp::base::DataVector alpha(gridStorage.getSize());

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < alpha.getSize(); i++) {
    //        alpha[i] = static_cast<double>(i);
    alpha[i] = dist(mt);
  }

  BOOST_TEST_MESSAGE("creating operation with unrefined grid");
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  sgpp::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());

  BOOST_TEST_MESSAGE("preparing operation for refined grid");
  eval->prepare();

  BOOST_TEST_MESSAGE("calculating result");

  double durationOverall = 0.0;
  size_t runs = 5;

  for (size_t i = 0; i < runs; i++) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    eval->mult(alpha, dataSizeVectorResult);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_TEST_MESSAGE("duration: " << elapsed_seconds.count());
    durationOverall += elapsed_seconds.count();
  }

  durationOverall /= static_cast<double>(runs);

  BOOST_TEST_MESSAGE("average duration: " << durationOverall);

  if (gridType == sgpp::base::GridType::Linear) {
    logger.outFile << datasetName << ","
                   << "Linear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    logger.outFile << datasetName << ","
                   << "ModLinear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else {
    throw nullptr;
  }
}

void prepareGrid(std::string fileName, sgpp::base::GridType gridType, size_t level) {
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = static_cast<int>(level);
  gridConfig.type_ = gridType;

  // Set Adaptivity
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 200;
  adaptivityConfig.numRefinements_ = 10;  // 6
  adaptivityConfig.percent_ = 100.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 50;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 50;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  std::string metaInformation =
      "refine: " + std::to_string(adaptivityConfig.numRefinements_) + " points: " +
      std::to_string(adaptivityConfig.numRefinementPoints_) + " iterations: " +
      std::to_string(SLESolverConfigRefine.maxIterations_);

  double lambda = 0.000001;

  bool verbose = true;
  sgpp::datadriven::MetaLearner* learner = new sgpp::datadriven::MetaLearner(
      gridConfig, SLESolverConfigRefine, SLESolverConfigFinal, adaptivityConfig, lambda, verbose);

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", true);
  parameters.addIDAttr("VERBOSE", true);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  //    parameters.addIDAttr("PLATFORM", "Intel(R) OpenCL");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  //    parameters.addIDAttr("WRITE_SOURCE", "true");
  //    parameters.addIDAttr("REUSE_SOURCE", "true");
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration;

  if (gridType == sgpp::base::GridType::Linear) {
    configuration = sgpp::datadriven::OperationMultipleEvalConfiguration(
        sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
        sgpp::datadriven::OperationMultipleEvalSubType::COMBINED);
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    configuration = sgpp::datadriven::OperationMultipleEvalConfiguration(
        sgpp::datadriven::OperationMultipleEvalType::STREAMING,
        sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP);
  }

  std::string content = uncompressFile(fileName);

  learner->learnString(configuration, content);

  BOOST_TEST_MESSAGE("info: grid preparation by metalearner is complete!");

  if (gridType == sgpp::base::GridType::Linear) {
    preparedGrids[fileName] = learner;
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    preparedGridsModLinear[fileName] = learner;
  } else {
    throw;
  }
}

void getRuntimeDataMining(sgpp::base::GridType gridType, const std::string& kernel,
                          std::string& fileName, std::string& datasetName, size_t level,
                          sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  if (gridType == sgpp::base::GridType::Linear) {
    if (preparedGrids.size() == 0) {
      prepareGrid(fileName, gridType, level);
    }
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    if (preparedGridsModLinear.size() == 0) {
      prepareGrid(fileName, gridType, level);
    }
  } else {
    throw;
  }

  sgpp::datadriven::MetaLearner* learner = nullptr;

  if (gridType == sgpp::base::GridType::Linear) {
    learner = preparedGrids[fileName];
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    learner = preparedGridsModLinear[fileName];
  } else {
    throw;
  }

  sgpp::base::Grid& grid = learner->getLearnedGrid();

  sgpp::base::GridStorage& gridStorage = grid.getStorage();

  BOOST_TEST_MESSAGE("size of refined grid: " << gridStorage.getSize());
  refinedGridSize = gridStorage.getSize();

  sgpp::base::DataVector alpha(gridStorage.getSize());

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < alpha.getSize(); i++) {
    //        alpha[i] = static_cast<double>(i);
    alpha[i] = dist(mt);
  }

  BOOST_TEST_MESSAGE("creating operation with unrefined grid");
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval(
      sgpp::op_factory::createOperationMultipleEval(grid, trainingData, configuration));

  sgpp::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());

  BOOST_TEST_MESSAGE("preparing operation for refined grid");
  eval->prepare();

  BOOST_TEST_MESSAGE("calculating result");

  double durationOverall = 0.0;
  size_t runs = 5;

  for (size_t i = 0; i < runs; i++) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    eval->mult(alpha, dataSizeVectorResult);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_TEST_MESSAGE("duration: " << elapsed_seconds.count());
    durationOverall += elapsed_seconds.count();
  }

  durationOverall /= static_cast<double>(runs);

  BOOST_TEST_MESSAGE("average duration: " << durationOverall);

  if (gridType == sgpp::base::GridType::Linear) {
    logger.outFile << datasetName << ","
                   << "Linear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    logger.outFile << datasetName << ","
                   << "ModLinear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else {
    throw nullptr;
  }
}

void getRuntimeDataMiningTransposed(
    sgpp::base::GridType gridType, const std::string& kernel, std::string& fileName,
    std::string& datasetName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  if (gridType == sgpp::base::GridType::Linear) {
    if (preparedGrids.size() == 0) {
      prepareGrid(fileName, gridType, level);
    }
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    if (preparedGridsModLinear.size() == 0) {
      prepareGrid(fileName, gridType, level);
    }
  } else {
    throw;
  }

  sgpp::datadriven::MetaLearner* learner = nullptr;

  if (gridType == sgpp::base::GridType::Linear) {
    learner = preparedGrids[fileName];
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    learner = preparedGridsModLinear[fileName];
  } else {
    throw;
  }

  sgpp::base::Grid& grid = learner->getLearnedGrid();

  sgpp::base::GridStorage& gridStorage = grid.getStorage();

  BOOST_TEST_MESSAGE("size of refined grid: " << gridStorage.getSize());
  refinedGridSize = gridStorage.getSize();

  sgpp::base::DataVector source(dataset.getNumberInstances());

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < source.getSize(); i++) {
    //        alpha[i] = static_cast<double>(i);
    source[i] = dist(mt);
  }

  BOOST_TEST_MESSAGE("creating operation with unrefined grid");
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval(
      sgpp::op_factory::createOperationMultipleEval(grid, trainingData, configuration));

  sgpp::base::DataVector gridSizeVectorResult(gridStorage.getSize());

  BOOST_TEST_MESSAGE("preparing operation for refined grid");
  eval->prepare();

  BOOST_TEST_MESSAGE("calculating result");

  double durationOverall = 0.0;
  size_t runs = 5;

  for (size_t i = 0; i < runs; i++) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    eval->multTranspose(source, gridSizeVectorResult);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_TEST_MESSAGE("duration: " << elapsed_seconds.count());
    durationOverall += elapsed_seconds.count();
  }

  durationOverall /= static_cast<double>(runs);

  BOOST_TEST_MESSAGE("average duration: " << durationOverall);

  if (gridType == sgpp::base::GridType::Linear) {
    logger.outFile << datasetName << ","
                   << "Linear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    logger.outFile << datasetName << ","
                   << "ModLinear"
                   << "," << kernel << "," << refinedGridSize << "," << durationOverall
                   << std::endl;
  } else {
    throw nullptr;
  }
}

void getRuntimeTransposed(sgpp::base::GridType gridType, const std::string& kernel,
                          std::string& fileName, std::string& datasetName, size_t level,
                          sgpp::base::AdaptivityConfiguration adaptivityConfig,
                          sgpp::datadriven::OperationMultipleEvalConfiguration configuration) {
  std::string content = uncompressFile(fileName);

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  sgpp::base::DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<sgpp::base::Grid> grid;

  if (gridType == sgpp::base::GridType::Linear) {
    grid = std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    grid = std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createModLinearGrid(dim));
  } else {
    throw;
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  BOOST_TEST_MESSAGE("dimensionality:        " << gridStorage.getDimension());

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(level);
  BOOST_TEST_MESSAGE("number of grid points: " << gridStorage.getSize());
  BOOST_TEST_MESSAGE("number of data points: " << dataset.getNumberInstances());

  doDirectedRefinements(adaptivityConfig, *grid, gridGen);
  //    doRandomRefinements(adaptivityConfig, *grid, gridGen);

  BOOST_TEST_MESSAGE("size of refined grid: " << gridStorage.getSize());
  refinedGridSize = gridStorage.getSize();

  sgpp::base::DataVector source(dataset.getNumberInstances());

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (size_t i = 0; i < source.getSize(); i++) {
    //        alpha[i] = static_cast<double>(i);
    source[i] = dist(mt);
  }

  BOOST_TEST_MESSAGE("creating operation with unrefined grid");
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval(
      sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  sgpp::base::DataVector gridSizeVectorResult(gridStorage.getSize());

  BOOST_TEST_MESSAGE("preparing operation for refined grid");
  eval->prepare();

  BOOST_TEST_MESSAGE("calculating result");

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  eval->multTranspose(source, gridSizeVectorResult);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  BOOST_TEST_MESSAGE("duration: " << elapsed_seconds.count());

  if (gridType == sgpp::base::GridType::Linear) {
    logger.outFile << datasetName << ","
                   << "Linear"
                   << "," << kernel << "," << refinedGridSize << "," << elapsed_seconds.count()
                   << std::endl;
  } else if (gridType == sgpp::base::GridType::ModLinear) {
    logger.outFile << datasetName << ","
                   << "ModLinear"
                   << "," << kernel << "," << refinedGridSize << "," << elapsed_seconds.count()
                   << std::endl;
  } else {
    throw nullptr;
  }
}

BOOST_AUTO_TEST_SUITE(HPCSE2015Linear)

BOOST_AUTO_TEST_CASE(StreamingDefault) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(1));
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
  parameters.addTextAttr("KERNEL_STORE_DATA", "array");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(1));
  parameters.addTextAttr("PLATFORM", "first");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementSteps[i];
    //        getRuntimeDataMining(GridType::Linear, "AVX", fileNames[i], datasetNames[i],
    //        levels[i], adaptivityConfig,
    //                configuration);
    getRuntimeDataMining(sgpp::base::GridType::Linear, "AVX", fileNames[i], datasetNames[i],
                         levels[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingSubspaceLinear) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
      sgpp::datadriven::OperationMultipleEvalSubType::COMBINED);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementSteps[i];
    //        getRuntimeDataMining(GridType::Linear, "Subspace", fileNames[i], datasetNames[i],
    //        levels[i], adaptivityConfig,
    //                configuration);
    getRuntimeDataMining(sgpp::base::GridType::Linear, "Subspace", fileNames[i], datasetNames[i],
                         levels[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingBase) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold = 0.0;

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::DEFAULT,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementSteps[i];
    //        getRuntimeDataMining(GridType::Linear, "Generic", fileNames[i], datasetNames[i],
    //        levels[i], adaptivityConfig,
    //                configuration);
    getRuntimeDataMining(sgpp::base::GridType::Linear, "Generic", fileNames[i], datasetNames[i],
                         levels[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingOCL) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(1));
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
  parameters.addTextAttr("KERNEL_STORE_DATA", "register");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementSteps[i];
    //        getRuntimeDataMining(GridType::Linear, "OCL (GPU)", fileNames[i], datasetNames[i],
    //        levels[i], adaptivityConfig,
    //                configuration);
    getRuntimeDataMining(sgpp::base::GridType::Linear, "OCL (GPU)", fileNames[i], datasetNames[i],
                         levels[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingOCLBlocking) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(4));
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(4));
  parameters.addTextAttr("KERNEL_STORE_DATA", "register");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementSteps[i];
    //        getRuntimeDataMining(GridType::Linear, "OCL blocked (GPU)", fileNames[i],
    //        datasetNames[i], levels[i],
    //                adaptivityConfig, configuration);
    getRuntimeDataMining(sgpp::base::GridType::Linear, "OCL blocked (GPU)", fileNames[i],
                         datasetNames[i], levels[i], configuration);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(HPCSE2015ModLinear)

BOOST_AUTO_TEST_CASE(StreamingBase) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::DEFAULT,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementStepsModLinear[i];
    //        getRuntimeDataMiningTransposed(GridType::ModLinear, "Generic", fileNames[i],
    //        datasetNames[i],
    //                levelsModLinear[i], adaptivityConfig, configuration);
    getRuntimeDataMiningTransposed(sgpp::base::GridType::ModLinear, "Generic", fileNames[i],
                                   datasetNames[i], levelsModLinear[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingOCL) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(1));
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
  parameters.addTextAttr("KERNEL_STORE_DATA", "register");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementStepsModLinear[i];
    //        getRuntimeDataMiningTransposed(GridType::ModLinear, "OCL (GPU)", fileNames[i],
    //        datasetNames[i],
    //                levelsModLinear[i], adaptivityConfig, configuration);
    getRuntimeDataMiningTransposed(sgpp::base::GridType::ModLinear, "OCL (GPU)", fileNames[i],
                                   datasetNames[i], levelsModLinear[i], configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingOCLFast) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(4));
  parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
  parameters.addIDAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", UINT64_C(8));
  parameters.addTextAttr("KERNEL_STORE_DATA", "register");
  parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCLFASTMP, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementStepsModLinear[i];
    //        getRuntimeDataMiningTransposed(GridType::ModLinear, "OCL blocked (GPU)", fileNames[i],
    //        datasetNames[i],
    //                levelsModLinear[i], adaptivityConfig, configuration);
    getRuntimeDataMiningTransposed(sgpp::base::GridType::ModLinear, "OCL blocked (GPU)",
                                   fileNames[i], datasetNames[i], levelsModLinear[i],
                                   configuration);
  }
}

BOOST_AUTO_TEST_CASE(StreamingOCLMask) {
  //    sgpp::base::AdaptivityConfiguration adaptivityConfig;
  //    adaptivityConfig.maxLevelType_ = false;
  //    adaptivityConfig.numRefinementPoints_ = REFINEMENT_POINTS;
  //    adaptivityConfig.percent_ = 200.0;
  //    adaptivityConfig.refinementThreshold_ = 0.0;

  sgpp::base::OCLOperationConfiguration parameters;
  parameters.addIDAttr("OCL_MANAGER_VERBOSE", false);
  parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.addTextAttr("PLATFORM", "NVIDIA CUDA");
  parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", UINT64_C(0));
  parameters.addIDAttr("MAX_DEVICES", UINT64_C(1));

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP, parameters);

  for (size_t i = 0; i < fileNames.size(); i++) {
    //        adaptivityConfig.numRefinements_ = refinementStepsModLinear[i];
    //        getRuntimeDataMiningTransposed(GridType::ModLinear, "OCL Mask (GPU)", fileNames[i],
    //        datasetNames[i],
    //                levelsModLinear[i], adaptivityConfig, configuration);
    getRuntimeDataMiningTransposed(sgpp::base::GridType::ModLinear, "OCL Mask (GPU)", fileNames[i],
                                   datasetNames[i], levelsModLinear[i], configuration);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
