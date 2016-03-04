// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
/*
 * DensityRegressionMetaLearner.hpp
 *
 *  Created on: Jan 7, 2016
 *      Author: pfandedd
 */

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

class PiecewiseConstantSmoothedRegressionMetaLearner {
 private:
  bool verbose;
  base::DataMatrix& dataset;
  base::DataVector& datasetValues;
  size_t dim;
  base::RegularGridConfiguration gridConfig;
  base::AdpativityConfiguration adaptConfig;
  solver::SLESolverConfiguration solverConfig;
  datadriven::RegularizationConfiguration regularizationConfig;

  /**
   * generates a regular grid
   * @param dim number of dimensions
   * @return grid
   */
  base::Grid* createRegularGrid(size_t dim);

  //    void splitset(std::vector<base::DataMatrix *>& strain, std::vector<base::DataVector *>&
  //  strainValues,
  // std::vector<base::DataMatrix *>& stest, std::vector<base::DataVector *>& stestValues,
  //  size_t kFold,
  // bool shuffleDataset, uint32_t shuffleSeed);

  void optimizeLambdaLog_(size_t kFold, size_t maxLevel,
                          double fastApproximationMSE,
                          size_t fastApproximationMaxLevel,
                          std::vector<base::DataMatrix>& trainingSets,
                          std::vector<base::DataVector>& trainingSetsValues,
                          std::vector<base::DataMatrix>& testSets,
                          std::vector<base::DataVector>& testSetsValues, size_t curLevel,
                          double lambdaLogStepSize,
                          double& bestLogLambda, double& bestMSE);

 public:
  PiecewiseConstantSmoothedRegressionMetaLearner(bool verbose,
      base::DataMatrix& trainingDataSet, base::DataVector& valuesDataSet,
      base::RegularGridConfiguration gridConfig,
      base::AdpativityConfiguration adaptConfig,
      solver::SLESolverConfiguration solverConfig,
      datadriven::RegularizationConfiguration regularizationConfig);

  double optimizeLambdaLog(size_t kFold, size_t maxLevel,
                            double fastApproximationMSE,
                            size_t fastApproximationMaxLevel);

  void optimizeLambdaLog(size_t kFold, size_t maxLevel,
                         double fastApproximationMSE,
                         size_t fastApproximationMaxLevel, std::shared_ptr<base::Grid>& bestGrid,
                         std::shared_ptr<base::DataVector>& bestAlpha, double& lambdaOpt);

  /**
   * Does the learning step on a given grid, training set and regularization parameter lambda
   *
   * @param train sample set
   * @param trainValues training values
   * @param lambda regularization parameter
   * @param fastApproximationMSE mse for stopping piecewise constant approximation tree creation
   * @param fastApproximationMaxLevel maximum level for stopping piecewise constant approximation
   *   tree creation
   * @param grid grid
   * @param alpha grid coefficients
   */
  void train(base::DataMatrix& train, base::DataVector& trainValues,
             double lambda, double fastApproximationMSE,
             size_t fastApproximationMaxLevel, std::shared_ptr<base::Grid>& grid,
             std::shared_ptr<base::DataVector>& alpha);

  double calculateMSE(base::Grid& grid, base::DataVector& alpha,
                       base::DataMatrix& testSubset,
                       base::DataVector& valuesTestSubset, bool verbose = false);
};

}  // namespace datadriven
}  // namespace sgpp

