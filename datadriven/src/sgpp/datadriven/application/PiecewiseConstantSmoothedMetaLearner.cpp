// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
/*
 * DensityRegressionMetaLearner.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: pfandedd
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/datadriven/application/PiecewiseConstantSmoothedMetaLearner.hpp>

#include <sgpp/datadriven/tools/DatasetTools.hpp>
#include <sgpp/datadriven/application/LearnerPiecewiseConstantSmoothedRegression.hpp>
#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/OperationPiecewiseConstantRegression.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

PiecewiseConstantSmoothedRegressionMetaLearner::PiecewiseConstantSmoothedRegressionMetaLearner(
  bool verbose, base::DataMatrix& trainingDataSet,
  base::DataVector& valuesDataSet, base::RegularGridConfiguration gridConfig,
  base::AdpativityConfiguration adaptConfig,
  solver::SLESolverConfiguration solverConfig,
  datadriven::RegularizationConfiguration regularizationConfig) :
  verbose(verbose), dataset(trainingDataSet), datasetValues(valuesDataSet),
  gridConfig(gridConfig), adaptConfig(
    adaptConfig), solverConfig(solverConfig),
  regularizationConfig(regularizationConfig) {
  this->dim = trainingDataSet.getNcols();
}

void PiecewiseConstantSmoothedRegressionMetaLearner::optimizeLambdaLog(
  size_t kFold, size_t maxLevel, double fastApproximationMSE,
  size_t fastApproximationMaxLevel, std::shared_ptr<base::Grid>& bestGrid,
  std::shared_ptr<base::DataVector>& bestAlpha, double& lambdaOpt) {
  lambdaOpt = this->optimizeLambdaLog(kFold, maxLevel, fastApproximationMSE,
                                      fastApproximationMaxLevel);
  this->train(dataset, datasetValues, lambdaOpt, fastApproximationMSE,
              fastApproximationMaxLevel, bestGrid, bestAlpha);
}

double PiecewiseConstantSmoothedRegressionMetaLearner::optimizeLambdaLog(
    size_t kFold, size_t maxLevel, double fastApproximationMSE,
    size_t fastApproximationMaxLevel) {
  std::vector<base::DataMatrix> trainingSets;
  std::vector<base::DataVector> trainingSetsValues;
  std::vector<base::DataMatrix> testSets;
  std::vector<base::DataVector> testSetsValues;

  DatasetTools::splitset(this->dataset, this->datasetValues, kFold, trainingSets,
                         trainingSetsValues, testSets,
                         testSetsValues);

  // initial values are pure dummy values
  double bestLambda = 0.0;
  double bestMSE = 0.0;

  this->optimizeLambdaLog_(kFold, maxLevel, fastApproximationMSE,
                           fastApproximationMaxLevel, trainingSets,
                           trainingSetsValues, testSets, testSetsValues, 0, 1.0,
                           bestLambda, bestMSE);
  bestLambda = pow(10.0, -bestLambda);

  if (verbose) {
    std::cout << "# -> bestLambda = " << bestLambda << std::endl;
    std::cout << "# -> bestMSE= " << bestMSE << std::endl;
  }

  return bestLambda;
}

void PiecewiseConstantSmoothedRegressionMetaLearner::optimizeLambdaLog_(size_t kFold,
    size_t maxLevel, double fastApproximationMSE,
    size_t fastApproximationMaxLevel, std::vector<base::DataMatrix>& trainingSets,
    std::vector<base::DataVector>& trainingSetsValues,
    std::vector<base::DataMatrix>& testSets,
    std::vector<base::DataVector>& testSetsValues, size_t curLevel,
    double lambdaLogStepSize,
    double& bestLogLambda, double& bestMSE) {
  if (verbose) {
    std::cout << "entering level=" << curLevel << " with lambda=" << pow(10.0,
              -bestLogLambda) << std::endl;
  }

  std::vector<double> logLambdaValues;

  if (curLevel == 0 && lambdaLogStepSize == 1.0) {
    for (size_t i = 4; i <= 10; i++) {
      logLambdaValues.push_back(static_cast<double>(i));
    }
  } else {
    logLambdaValues.push_back(bestLogLambda - lambdaLogStepSize);

    if (curLevel == 0) {
      logLambdaValues.push_back(bestLogLambda);
    }

    logLambdaValues.push_back(bestLogLambda + lambdaLogStepSize);
  }

  bool firstValue = true;

  for (double curLogLambda : logLambdaValues) {
    double curLambda = pow(10, -curLogLambda);
    std::cout << "curLambda: " << curLambda << std::endl;

    // cross-validation
    double curMeanMSE = 0.0;

    for (size_t j = 0; j < kFold; j++) {
      std::shared_ptr<base::Grid> grid;
      std::shared_ptr<base::DataVector> alpha;
      // compute density
      train(trainingSets[j], trainingSetsValues[j], curLambda, fastApproximationMSE,
            fastApproximationMaxLevel,
            grid, alpha);

      double mse = this->calculateMSE(*grid, *alpha, testSets[j], testSetsValues[j]);
      curMeanMSE += mse;
    }

    curMeanMSE /= static_cast<double>(kFold);

    if ((curLevel == 0 && firstValue) || curMeanMSE < bestMSE) {
      bestMSE = curMeanMSE;
      bestLogLambda = curLogLambda;
      firstValue = false;

      if (verbose) {
        std::cout << "new best lambda!" << std::endl;
      }
    } else {
      if (curLevel == 0) {
        break;
      }
    }

    if (verbose) {
      std::cout << "# lambda: " << curLambda << " curMeanMSE: " << curMeanMSE <<
                " bestLambda: "
                << pow(10.0, -bestLogLambda) << " bestMSE: " << bestMSE <<
                " lambdaLogStepSize: "
                << lambdaLogStepSize << std::endl;
    }
  }

  if (curLevel < maxLevel) {
    this->optimizeLambdaLog_(kFold, maxLevel, fastApproximationMSE,
                             fastApproximationMaxLevel, trainingSets,
                             trainingSetsValues, testSets, testSetsValues, curLevel + 1,
                             lambdaLogStepSize / 2.0, bestLogLambda,
                             bestMSE);
  }
}

void PiecewiseConstantSmoothedRegressionMetaLearner::train(
    base::DataMatrix& train, base::DataVector& trainValues, double lambda,
    double fastApproximationMSE, size_t fastApproximationMaxLevel,
    std::shared_ptr<base::Grid>& grid,
    std::shared_ptr<base::DataVector>& alpha) {
  sgpp::datadriven::OperationPiecewiseConstantRegression
  piecewiseRegressorOperator(train, trainValues);

  if (verbose) {
    std::cout << "creating piecewise-constant representation..." << std::endl;
  }

  //    std::unique_ptr<sgpp::datadriven::HistogramTree::Node> piecewiseRegressor =
  //   piecewiseRegressorOperator.hierarchize(
  //            0.0, 30);
  std::unique_ptr<sgpp::datadriven::PiecewiseConstantRegression::Node>
  piecewiseRegressor = piecewiseRegressorOperator.hierarchize(
                         fastApproximationMSE, fastApproximationMaxLevel);

  if (verbose) {
    std::cout << "piecewise-constant representation mse: " <<
              piecewiseRegressor->getMSE() << std::endl;
  }

  if (verbose) {
    std::cout << "piecewise-constant representation created, smoothing..." <<
              std::endl;
  }

  sgpp::datadriven::LearnerPiecewiseConstantSmoothedRegression learner(gridConfig,
      adaptConfig, solverConfig, regularizationConfig,
      true);

  // initialize standard grid and alpha vector
  grid = std::shared_ptr<base::Grid>(createRegularGrid(this->dim));
  alpha = std::make_shared<base::DataVector>(grid->getSize());

  learner.train(*piecewiseRegressor, *grid, *alpha, lambda);

  if (verbose) {
    std::cout << "smoothing finished" << std::endl;
  }
}

base::Grid* PiecewiseConstantSmoothedRegressionMetaLearner::createRegularGrid(
  size_t dim) {
  base::Grid* grid = nullptr;

  // load grid
  if (gridConfig.type_ == base::GridType::Linear) {
    grid = base::Grid::createLinearGrid(dim);
  } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
    grid = base::Grid::createLinearBoundaryGrid(dim, 0);
  } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
    grid = base::Grid::createLinearBoundaryGrid(dim);
  } else {
    throw base::application_exception("DensityRegressionMetaLearner::initialize : grid type is "
      "not supported");
  }

  grid->getGenerator().regular(gridConfig.level_);

  return grid;
}

double PiecewiseConstantSmoothedRegressionMetaLearner::calculateMSE(
  base::Grid& grid, base::DataVector& alpha,
  base::DataMatrix& testSubset, base::DataVector& valuesTestSubset,
  bool verbose) {
  double mse = 0.0;

  std::unique_ptr<base::OperationEval> opEval(sgpp::op_factory::createOperationEval(grid));

  for (size_t i = 0; i < testSubset.getNrows(); i++) {
    base::DataVector point(dim);
    testSubset.getRow(i, point);
    double approximation = opEval->eval(alpha, point);
    mse += (approximation - valuesTestSubset[i]) * (approximation -
           valuesTestSubset[i]);

    if (verbose && i < 100) {
      std::cout << "mine: " << approximation << " reference: " << valuesTestSubset[i]
                << std::endl;
    }
  }

  return mse;
}

}  // namespace datadriven
}  // namespace sgpp

