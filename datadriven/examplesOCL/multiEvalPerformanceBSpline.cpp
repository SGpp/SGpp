// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <random>
#include <string>

#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/operation/hash/common/basis/BsplineBasis.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"

SGPP::float_t calculateMSE(SGPP::base::DataVector vector1,
                           SGPP::base::DataVector vector2) {
  double mse = 0.0;

  for (size_t i = 0; i < vector1.getSize(); i++) {
    mse += (vector1[i] - vector2[i]) * (vector1[i] - vector2[i]);
  }

  mse /= static_cast<double>(vector1.getSize());
  return mse;
}

int main(int argc, char** argv) {
  std::string fileName = "friedman_4d.arff";
  //  std::string fileName = "debugging_small.arff";

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& trainingData = dataset.getData();

  // create a two-dimensional piecewise bi-linear grid
  const size_t dim = dataset.getDimension();
  const size_t degree = 5;
  SGPP::base::Grid* grid = SGPP::base::Grid::createBsplineGrid(dim, degree);
  SGPP::base::GridStorage* gridStorage = grid->getStorage();
  std::cout << "dimensionality:        " << gridStorage->dim() << "\n";
  // create regular grid, level 3
  const SGPP::base::GridIndex::level_type level = 8;
  SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  std::cout << "number of grid points: " << gridStorage->size() << "\n";

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1.0, 100.0);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL);

  SGPP::base::OperationMultipleEval* eval =
      SGPP::op_factory::createOperationMultipleEval(*grid, trainingData,
                                                    configuration);

  {
    // create coefficient vector
    SGPP::base::DataVector alpha(gridStorage->size());

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alpha[i] = dist(mt);
    }

    std::cout << "OCL B-spline evaluation...\n";
    SGPP::base::DataVector result(dataset.getNumberInstances());
    eval->eval(alpha, result);

    /*std::cout << "calculating comparison values...\n";
     SGPP::base::OperationNaiveEval* evalCompare =
     SGPP::op_factory::createOperationNaiveEval(*grid);
     SGPP::base::DataVector evaluationPoint(dim);
     SGPP::base::DataVector resultCompare(dataset.getNumberInstances());

     for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
     trainingData->getRow(i, evaluationPoint);
     resultCompare[i] = evalCompare->eval(alpha, evaluationPoint);
     }

     std::cout << "mse: " << calculateMSE(result, resultCompare) << "\n\n";*/
  }

  {
    // create coefficient vector
    SGPP::base::DataVector alpha(dataset.getNumberInstances());

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alpha[i] = dist(mt);
    }

    std::cout << "OCL B-spline transpose multiplication...\n";
    SGPP::base::DataVector result(gridStorage->size());
    eval->multTranspose(alpha, result);

    /*std::cout << "calculating comparison values...\n";
     SGPP::base::SBsplineBase base(degree);
     SGPP::base::DataVector evaluationPoint(dim);
     SGPP::base::DataVector resultCompare(gridStorage->size());

     for (size_t j = 0; j < gridStorage->size(); j++) {
     const SGPP::base::GridIndex& gp = *gridStorage->get(j);
     resultCompare[j] = 0.0;

     for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
     SGPP::float_t curVal = 1.0;
     trainingData->getRow(i, evaluationPoint);

     for (size_t t = 0; t < dim; t++) {
     SGPP::float_t val1d = base.eval(gp.getLevel(t), gp.getIndex(t),
     evaluationPoint[t]);

     if (val1d == 0.0) {
     curVal = 0.0;
     break;
     }

     curVal *= val1d;
     }

     resultCompare[j] += alpha[i] * curVal;
     }
     }

     std::cout << "mse: " << calculateMSE(result, resultCompare) << "\n";*/
  }
}
