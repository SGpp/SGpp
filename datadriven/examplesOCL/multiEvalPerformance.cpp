// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of
// distribution and use, please see the copyright notice
// provided with SG++ or at sgpp.sparsegrids.org.

#include <random>
#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"

int main(int argc, char** argv) {
  //  std::string fileName = "friedman_4d_2000.arff";
  std::string fileName = "debugging.arff";

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  // SGPP::base::DataVector *classes = dataset.getClasses();
  SGPP::base::DataMatrix& trainingData = dataset.getData();

  // create a two-dimensional piecewise bi-linear grid
  size_t dim = dataset.getDimension();
  std::unique_ptr<SGPP::base::Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  SGPP::base::GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;
  // create regular grid, level 3
  uint32_t level = 4;
  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage->size() << std::endl;

  // create coefficient vector
  SGPP::base::DataVector alpha(gridStorage->size());
  alpha.setAll(0.0);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  for (unsigned int i = 0; i < alpha.getSize(); i++) {
    alpha[i] = dist(mt);
    //    std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;
  }

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL);

  SGPP::base::OperationMultipleEval* eval =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData,
                                                    configuration);

  SGPP::base::DataVector result(dataset.getNumberInstances());

  eval->eval(alpha, result);

  //  std::cout << "result: ";
  //  for (size_t i = 0; i < result.getSize(); i++) {
  //    if (i > 0)
  //      std::cout << ",";
  //    std::cout << result[i];
  //  }
  //  std::cout << std::endl;

  std::cout << "calculating comparison values..." << std::endl;

  SGPP::base::OperationMultipleEval* evalCompare =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);
  SGPP::base::DataVector resultCompare(dataset.getNumberInstances());
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
