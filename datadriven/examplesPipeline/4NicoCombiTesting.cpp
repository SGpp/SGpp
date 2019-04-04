/*
 * 4NicoCombiTesting.cpp
 *
 *  Created on: Feb 14, 2019
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
  size_t dim, level, refinements;

  if (argc > 2) {
    dim = std::stoi(argv[1]);
    level = std::stoi(argv[2]);
  } else {
    cout << "setting defaults: dim = 2, level = 3" << std::endl;
    dim = 2;
    level = 3;
  }
  if (argc > 3) {
    refinements = std::stoi(argv[3]);
  } else {
    cout << "setting default number for refinements: 10" << std::endl;
    refinements = 10;
  }

  sgpp::datadriven::CombiConfigurator configurator;
  configurator.initAdaptiveScheme(dim, level);

  vector<combiConfig> vec;

  for (size_t a = 1; a <= refinements; a++) {
    cout << "REFINEMENT " << a << std::endl;
    configurator.getCombiScheme(vec);
    size_t i = 0;
    while (i + 1 < vec.size() && !configurator.isRefinable(vec.at(i))) {
      i++;
    }
    configurator.refineComponent(vec.at(i));
    configurator.getCombiScheme(vec);
  }

  /*
    for (int i = 0; i < vec.size(); i++) {
      if (configurator.isRefinable(vec.at(i))) {
        cout << "TRUE" << std::endl;
      } else {
        cout << "FALSE" << std::endl;
      }
    }
  */
  cout << "###################################################" << std::endl
       << "TEST \"4NicoCombiTesting\" TERMINATED SUCCESSFULLY!" << std::endl
       << "###################################################" << std::endl;
}
