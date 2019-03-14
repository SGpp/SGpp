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
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <string>
#include <vector>

using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using std::cout;
using std::string;
using std::to_string;
using std::vector;

string getTime() {
  time_t now = time(0);
  std::tm *time = std::localtime(&now);
  return to_string(1900 + time->tm_year) + "_" + to_string(1 + time->tm_mon) + "_" +
         to_string(time->tm_mday) + "_" + to_string(time->tm_hour) + "_" + to_string(time->tm_min) +
         "_" + to_string(time->tm_sec);
}

int main(int argc, char **argv) {
  cout << "Starting the evaluation. Time: " << getTime() << std::endl;
  std::ofstream myfile;
  string titel = "evaluation" + getTime();
  myfile.open(titel);
  myfile << "EVALUATION " << getTime() << std::endl;

  vector<string> paths;
  paths.push_back("classRipley");
  paths.push_back("classIris");
  paths.push_back("classDR10");

  ClassificationMinerFactory factory;
  for (string path : paths) {
    freopen((titel + "_" + path).c_str(), "w", stdout);
    cout << path;
    myfile << std::endl << std::endl << getTime() << path << std::endl;
    auto miner = std::unique_ptr<SparseGridMiner>(
        factory.buildMiner("classificationConfigs/" + path + ".json"));
    miner->learn(true);
    std::cout << std::endl;
  }

  myfile.close();
  return 0;
}
