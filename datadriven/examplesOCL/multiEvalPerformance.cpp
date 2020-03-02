// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

#include <random>
#include <string>


int main(int argc, char** argv) {
  //  std::string fileName = "friedman_4d_2000.arff";
  std::string fileName = "debugging.arff";

  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromFile(fileName);

  // sgpp::base::DataVector *classes = dataset.getClasses();
  sgpp::base::DataMatrix& trainingData = dataset.getData();

  // create a two-dimensional piecewise bi-linear grid
  size_t dim = dataset.getDimension();
  std::unique_ptr<sgpp::base::Grid> grid =
      std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;
  // create regular grid, level 3
  uint32_t level = 4;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;

  // create coefficient vector
  sgpp::base::DataVector alpha(gridStorage.getSize());

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (unsigned int i = 0; i < alpha.getSize(); i++) {
    alpha[i] = dist(mt);
    //    std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;
  }

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCL);

  std::unique_ptr<sgpp::base::OperationMultipleEval> eval =
      std::unique_ptr<sgpp::base::OperationMultipleEval>(
          sgpp::op_factory::createOperationMultipleEval(*grid, trainingData, configuration));

  sgpp::base::DataVector result(dataset.getNumberInstances());

  eval->eval(alpha, result);

  //  std::cout << "result: ";
  //  for (size_t i = 0; i < result.getSize(); i++) {
  //    if (i > 0)
  //      std::cout << ",";
  //    std::cout << result[i];
  //  }
  //  std::cout << std::endl;

  std::cout << "calculating comparison values..." << std::endl;

  std::unique_ptr<sgpp::base::OperationMultipleEval> evalCompare =
      std::unique_ptr<sgpp::base::OperationMultipleEval>(
          sgpp::op_factory::createOperationMultipleEval(*grid, trainingData));
  sgpp::base::DataVector resultCompare(dataset.getNumberInstances());
  evalCompare->eval(alpha, resultCompare);

  //  std::cout << "resultCompare: ";
  //  for (size_t i = 0; i < resultCompare.getSize(); i++) {
  //    if (i > 0)
  //      std::cout << ",";
  //    std::cout << resultCompare[i];
  //  }
  //  std::cout << std::endl;

  //  std::cout << "diff: ";
  //  for (size_t i = 0; i < result.getSize(); i++) {
  //    if (i > 0)
  //      std::cout << ",";
  //    std::cout << (result[i] - resultCompare[i]);
  //  }
  //  std::cout << std::endl;

  double mse = 0.0;

  for (size_t i = 0; i < result.getSize(); i++) {
    std::cout << "comp: " << (result[i] - resultCompare[i]) << std::endl;
    mse += (result[i] - resultCompare[i]) * (result[i] - resultCompare[i]);
  }

  mse = mse / static_cast<double>(result.getSize());
  std::cout << "mse: " << mse << std::endl;
}
