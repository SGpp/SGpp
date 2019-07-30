/* Copyright (C) 2008-today The SG++ project
 *
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 */

#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::DataSourceBuilder;
using sgpp::datadriven::DataSourceCrossValidation;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::ModelFittingLeastSquares;
using sgpp::datadriven::FitterConfigurationLeastSquares;
using sgpp::datadriven::MSE;
using sgpp::datadriven::Scorer;
using sgpp::datadriven::SparseGridMinerCrossValidation;
using sgpp::base::GridType;

int main(int argc, char** argv) {
  /**
   * use immediately invoked lambda expression to get the path to a configuration file.
   */
  const std::string path = [argc, &argv]() {
    if (argc != 2) {
      std::cout << "No or bad path given, aborting\n";
      exit(1);
      return std::string{};
    } else {
      return std::string{argv[1]};
    }
  }();

  /**
   * In order to read a dataset from disk, we need an instance of a #sgpp::datadriven::DataSource
   * object, that is constructed using a builder pattern. Since we only want to read a from disk and
   * use all samples it provides, we only pass in the path. Everything else is managed by default
   * values and auto detection of extensions.
   */
  DataSourceCrossValidation *dataSource =
      DataSourceBuilder().withPath(path).crossValidationAssemble();
  std::cout << "reading input file: " << path << std::endl;

  /**
   * We want to perform least squares regression now. First we need to set up our Fitter using a
   * configuration structure.
   */
  auto config = FitterConfigurationLeastSquares{};

 /**
   * We first set up the provided default parameters enabled for least squares regression
   */
  config.setupDefaults();
  /**
   * Everything that does not match the default values is then adapted.
   */
  auto& gridConfig = config.getGridConfig();
  gridConfig.level_ = 2;
  gridConfig.type_ = GridType::ModLinear;
  gridConfig.dim_ = dataSource->getNextSamples()->getDimension();
  auto& regularizationConfig = config.getRegularizationConfig();
  regularizationConfig.lambda_ = 10e-1;

  /**
   * Based on our configuration, we then can create a fitter object.
   */
  ModelFittingLeastSquares *model = new ModelFittingLeastSquares(config);
  /**
   * We want to perform 5 Fold cross validation on our model. To assess the quality of the
   * regression algorithm, we use the mean squared error (MSE) as an error metric. To ensure testing
   * and training data are not taken from an ordered distribution, we will permute the values that
   * go into testing and training dataset.
   */
  Scorer *scorer = new Scorer(new MSE{});

  /**
   * Create a sparse grid miner that performs cross validation. The number of folds is 5 per
   * default.
   */
  SparseGridMinerCrossValidation miner(dataSource, model, scorer);
  /**
   * Here the actual learning process is launched. The miner will perform k-fold cross validation
   * and print the mean score as well as the standard deviation.
   */
  miner.learn(true);
  return 0;
}
