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
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SimpleSplittingScorer.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <memory>
#include <string>

using sgpp::datadriven::FileBasedDataSourceBuilder;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::ModelFittingLeastSquares;
using sgpp::datadriven::DataMiningConfigurationLeastSquares;
using sgpp::datadriven::Metric;
using sgpp::datadriven::MSE;
using sgpp::datadriven::SimpleSplittingScorer;

int main(int argc, char **argv) {
  // input
  std::string path;
  if (argc != 2) {
    std::cout << "No or bad path given, aborting" << std::endl;
    exit(-1);
  } else {
    path = std::string(argv[1]);
  }

  auto dataSource = FileBasedDataSourceBuilder().withPath(path).assemble();
  std::cout << "reading input file: " << path << std::endl;
  auto dataset = dataSource->getNextSamples();

  double threshold = 0;
  size_t maxRefinenum = 4;

  // regression
  auto config = std::make_shared<DataMiningConfigurationLeastSquares>();
  config->addIDAttr("trainPortion", 0.5);
  config->addIDAttr("seed", static_cast<uint64_t>(42));

  auto gridConfig = config->getGridConfig();
  gridConfig.dim_ = dataset->getDimension();
  config->setGridConfig(gridConfig);

  std::cout << "fitting model" << std::endl;
  auto fitter = std::make_shared<ModelFittingLeastSquares>(config);
  fitter->fit(*dataset);

  // 4. load simple splitting scorer
  std::shared_ptr<Metric> metric = std::make_shared<MSE>();
  SimpleSplittingScorer scorer(metric, fitter, *config);

  double score = scorer.getScore(*dataset);
  std::cout << "MSE is " << score << std::endl;
  size_t iter = 0;
  while (score > threshold && iter < maxRefinenum) {
    std::cout << "refining grid" << std::endl;
    fitter->refine();
    std::cout << "updating grid" << std::endl;
    fitter->update(*dataset);
    score = scorer.getScore(*dataset);
    std::cout << "MSE of refined grid is " << score << std::endl;
    iter++;
  }
  return 0;
}
