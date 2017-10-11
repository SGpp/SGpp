/* Copyright (C) 2008-today The SG++ project

 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * regressionPipeline.cpp
 *
 *  Created on: 01.06.2016
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::DataSourceBuilder;
using sgpp::datadriven::DataSource;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::ModelFittingLeastSquares;
using sgpp::datadriven::FitterConfigurationLeastSquares;
using sgpp::datadriven::MSE;
using sgpp::datadriven::Scorer;
using sgpp::datadriven::CrossValidation;
using sgpp::datadriven::SplittingScorer;
using sgpp::datadriven::RandomShufflingFunctor;
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
  auto dataSource = std::unique_ptr<DataSource>(DataSourceBuilder().withPath(path).assemble());
  std::cout << "reading input file: " << path << std::endl;
  /**
   * Once we have a data source, we can read the contents of the stored dataset.
   */
  auto dataset = std::unique_ptr<Dataset>(dataSource->getNextSamples());

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
  gridConfig.dim_ = dataset->getDimension();
  auto& regularizationConfig = config.getRegularizationConfig();
  regularizationConfig.lambda_ = 10e-1;

  std::cout << "starting 5 fold cross validation with seed 42" << std::endl;
  /**
   * Based on our configuration, we then can create a fitter object.
   */
  auto model = ModelFittingLeastSquares{config};
  double stdDeviation = 0;
  /**
   * We want to perform 5 Fold cross validation on our model. To assess the quality of the
   * regression algorithm, we use the mean squared error (MSE) as an error metric. To ensure testing
   * and training data are not taken from an ordered distribution, we will permute the values that
   * go into testing and training dataset and use a random seed of 42 for shuffling.
   */
  CrossValidation scorer{new MSE{}, new RandomShufflingFunctor{}, 42, 5};

  /**
   * Here the actual learning process is launched. The calculate Score method will perform 5 fold
   * cross validation using our least squares model and return the average MSE and the standard
   * deviation of the individual fold MSEs.
   */
  auto score = scorer.calculateScore(model, *dataset, &stdDeviation);

  std::cout << "Score = " << score << " with stdDeviation " << stdDeviation << std::endl;

  return 0;
}
