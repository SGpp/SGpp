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

#include <sgpp/datadriven/datamining/builder/FileBasedDataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::FileBasedDataSourceBuilder;
using sgpp::datadriven::Dataset;

using sgpp::datadriven::ModelFittingLeastSquares;
using sgpp::datadriven::DataMiningConfigurationLeastSquares;

int main(int argc, char **argv) {
  // input
  auto dataSource =
      FileBasedDataSourceBuilder()
          .withPath(
              "/home/michael/Projects/SGpp/datasets/bupa_liver/liver-disorders_normalized.arff.gz")
          .assemble();
  std::cout << "reading input file: "
            << "liver-disorders_normalized.arff.gz" << std::endl;
  auto dataset = dataSource->getNextSamples();

  // regression
  DataMiningConfigurationLeastSquares config;

  auto gridConfig = config.getGridConfig();
  gridConfig.dim_ = dataset->getDimension();
  config.setGridConfig(gridConfig);

  auto fitter = std::make_unique<ModelFittingLeastSquares>(config);
  fitter->fit(*(dataset.get()));
  std::cout << "fitting model" << std::endl;

  auto surplusses = fitter->getSurpluses();

  std::cout << "surplusses are: " << surplusses->toString();

  return 0;
}
