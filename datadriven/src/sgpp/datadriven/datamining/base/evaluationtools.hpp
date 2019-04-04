/*
 * evaluationtools.hpp
 *
 *  Created on: Mar 14, 2019
 *      Author: nico
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <ctime>
#include <string>
#include <vector>
using std::string;
using std::to_string;
using std::vector;
namespace evalu {
string getTime();

string getMSek();

string dmtS(sgpp::base::DataMatrix &m);

string mtS(vector<vector<size_t>> &m);
}
