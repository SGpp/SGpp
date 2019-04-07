/*
 * 0_Evaluation.cpp
 *
 *  Created on: Mar 13, 2019
 *      Author: nico
 */

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/base/evaluationtools.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <sstream>
#include <string>
#include <vector>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using std::cout;
using std::string;
using std::to_string;
using std::vector;
using evalu::getTime;

int main(int argc, char **argv) {
  cout << "Starting the evaluation. Time: " << getTime() << std::endl;
  cout << evalu::getMSek();
  std::ofstream myfile;
  string titel = "evaluation" + getTime();
  myfile.open(titel);
  myfile << "EVALUATION " << getTime() << std::endl;

  vector<string> paths;
  for (size_t i = 0; i < 25; i++) {
    paths.push_back("e1_Circles500_component");
    paths.push_back("e1_Circles500_regular");
  }
  for (size_t i = 0; i < 0; i++) {
    paths.push_back("classIris_component");
    paths.push_back("classIris_regular");
  }
  for (size_t i = 0; i < 0; i++) {
    paths.push_back("classDR10_component");
    paths.push_back("classDR10_regular");
  }

  ClassificationMinerFactory factory;
  for (string path : paths) {
    // timestamp in main evaluation file
    myfile << std::endl << std::endl << getTime() << " " << evalu::getMSek() << path << std::endl;

    // writing in the child file
    freopen((titel + "_" + path + evalu::getTime() + evalu::getMSek()).c_str(), "w", stdout);
    cout << path << std::endl << std::endl;
    // copying the config.json
    std::ifstream inFile;
    inFile.open("classificationConfigs/" + path + ".json");
    string line;
    while (std::getline(inFile, line)) {
      cout << line;
    }
    cout << std::endl << std::endl;

    auto miner = std::unique_ptr<SparseGridMiner>(
        factory.buildMiner("classificationConfigs/" + path + ".json"));
    miner->learn(true);
    std::cout << std::endl << "THIS EVALUATION RUN TERMINATED JUST AS PLANED";
  }
  myfile << std::endl << "DONE" << std::endl << getTime() << std::endl << evalu::getMSek();
  myfile.close();
  return 0;
}
