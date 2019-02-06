/*
 * CombiTesting.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: nico
 */

#include <iostream>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombiGrid.hpp>

using std::cout;
using std::vector;
using sgpp::datadriven::combiConfig;

int main(int argc, char** argv) {
  const std::string path = [argc, &argv]() {
    if (argc != 2) {
      std::cout << "No or bad path given, aborting\n";
      exit(1);
      return std::string{};
    } else {
      return std::string{argv[1]};
    }
  }();
  sgpp::datadriven::CombiConfigurator combi;
  vector<sgpp::datadriven::combiConfig> combiconfig;
  combi.getStandardCombi(combiconfig, 2, 5);

  const sgpp::datadriven::DataMiningConfigParser parser(path);

  sgpp::datadriven::FitterConfigurationDensityEstimation config{};
  config.readParams(parser);
  auto fittingbase = new sgpp::datadriven::ModelFittingDensityEstimationCombiGrid(config);
  cout << "ModelFittingDensityEstimationCombiGrid created!" << std::endl;
  sgpp::base::DataMatrix data;
  fittingbase->fit(data);
  return 0;
}
