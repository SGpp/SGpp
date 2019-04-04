/*
 * 2NicoCombiTesting.cpp
 *
 *  Created on: Feb 13, 2019
 *      Author: nico
 */

#include <iostream>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <stdexcept>
#include <string>

using std::cout;
using std::vector;
using sgpp::datadriven::combiConfig;

int main(int argc, char** argv) {
  size_t dim, level;

  if (argc > 2) {
    dim = std::stoi(argv[1]);
    level = std::stoi(argv[2]);
  } else {
    cout << "setting defaults: dim = 2, level = 3" << std::endl;
    dim = 2;
    level = 3;
  }

  sgpp::datadriven::CombiConfigurator configurator;
  configurator.initAdaptiveScheme(dim, level);

  vector<combiConfig> vec;

  configurator.getCombiScheme(vec);

  for (int i = 0; i < vec.size(); i++) {
    if (configurator.isRefinable(vec.at(i))) {
      cout << "TRUE" << std::endl;
    } else {
      cout << "FALSE" << std::endl;
    }
  }
  cout << std::endl << "Test \"2NicoCombiTesting\" terminated successfully!" << std::endl;
}
