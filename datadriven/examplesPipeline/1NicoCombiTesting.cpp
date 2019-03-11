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
#include "../src/sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp"

using std::cout;
using std::vector;
using sgpp::datadriven::combiConfig;

int main(int argc, char** argv) {
  cout << "Standard combi for dim = 3, level = 3\n";
  vector<sgpp::datadriven::combiConfig> combiconfig;
  sgpp::datadriven::CombiConfigurator combiconfigurator;
  combiconfigurator.getStandardCombi(combiconfig, 3, 3);

  cout << "Printing CombiConfig: \n";
  for (auto v : combiconfig) {
    cout << v.coef << "[";
    for (auto b : v.levels) {
      cout << b << " ";
    }
    cout << "]\n";
  }

  combiconfig.clear();
  sgpp::datadriven::CombiConfigurator configurator;
  configurator.test(combiconfig);
  cout << "Printing CombiConfig: \n";
  for (auto v : combiconfig) {
    cout << v.coef << "[";
    for (auto b : v.levels) {
      cout << b << " ";
    }
    cout << "]\n";
  }
  return 0;
}
