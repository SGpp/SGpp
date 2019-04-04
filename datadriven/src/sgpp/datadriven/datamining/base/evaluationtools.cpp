/*
 * evaluationtools.cpp
 *
 *  Created on: Mar 14, 2019
 *      Author: nico
 */

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/datamining/base/evaluationtools.hpp>

#include <chrono>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::to_string;
using std::vector;
namespace evalu {
string getTime() {
  time_t now = time(0);
  std::tm *time = std::localtime(&now);
  return to_string(1900 + time->tm_year) + "_" + to_string(1 + time->tm_mon) + "_" +
         to_string(time->tm_mday) + "_" + to_string(time->tm_hour) + "_" + to_string(time->tm_min) +
         "_" + to_string(time->tm_sec);
}

string getMSek() {
  size_t sek = std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::system_clock::now().time_since_epoch())
                   .count();
  return "MSECS: " + to_string(sek) + "\n";
}

string dmtS(sgpp::base::DataMatrix &m) {
  std::stringstream s;
  for (size_t i = 0; i < m.getNcols(); i++) {
    s << "[";
    for (size_t j = 0; j < m.getNrows(); j++) {
      s << static_cast<int>(m.get(j, i));
      s << ", ";
    }
    s << "]\n";
  }
  return s.str();
}

string mtS(vector<vector<size_t>> &m) {
  string s = "";
  for (size_t i = 0; i < m.size(); i++) {
    vector<size_t> v;
    s += "[";
    for (size_t j = 0; j < m.at(i).size(); j++) {
      s += m.at(i).at(j) + " ,";
    }
    s += "]\n";
  }
  return s;
}
}
