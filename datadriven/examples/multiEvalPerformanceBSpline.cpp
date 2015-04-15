/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <random>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

int main(int argc, char** argv) {

  std::string fileName = "friedman_4d.arff";
  //  std::string fileName = "debugging_small.arff";

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

  // create a two-dimensional piecewise bi-linear grid
  const size_t dim = dataset.getDimension();
  const size_t degree = 3;
  SGPP::base::Grid* grid = SGPP::base::Grid::createBsplineGrid(dim, degree);
  SGPP::base::GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;
  // create regular grid, level 3
  const SGPP::base::GridIndex::level_type level = 6;
  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage->size() << std::endl;

  // create coefficient vector
  SGPP::base::DataVector alpha(gridStorage->size());
  alpha.setAll(0.0);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1.0, 100.0);

  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha[i] = dist(mt);
    //    std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;
  }

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
  configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
  configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

  SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

  //  SGPP::base::OperationMultipleEval *eval =
  //  SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

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

  SGPP::base::OperationNaiveEval* evalCompare =
    SGPP::op_factory::createOperationNaiveEval(*grid);
  SGPP::base::DataVector evaluationPoint(dim);
  SGPP::base::DataVector resultCompare(dataset.getNumberInstances());

  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    trainingData->getRow(i, evaluationPoint);
    resultCompare[i] = evalCompare->eval(alpha, evaluationPoint);
  }

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
    //std::cout << "comp: " << (result[i] - resultCompare[i]) << std::endl;
    mse += (result[i] - resultCompare[i]) * (result[i] - resultCompare[i]);
  }

  mse /= static_cast<double>(result.getSize());
  std::cout << "mse: " << mse << std::endl;
}

