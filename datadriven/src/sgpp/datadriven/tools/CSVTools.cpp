// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
// Created on: 18.12.2017
// Author: Eric Koepke

#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <sgpp/base/exception/file_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>

namespace sgpp {
namespace datadriven {

Dataset CSVTools::readCSV(const std::string& filename, bool skipFirstLine) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  size_t numberInstances;
  size_t dimension;
  size_t instanceNo = 0;

  readCSVSize(filename, numberInstances, dimension);
  if (skipFirstLine && !myfile.eof()) {
    numberInstances--;
    std::getline(myfile, line);
  }
  Dataset dataset(numberInstances, dimension);

  while (!myfile.eof()) {
    std::getline(myfile, line);

    if (!line.empty()) {
      writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }
  }

  myfile.close();

  return dataset;
}

void CSVTools::readCSVSize(const std::string& filename,
                             size_t& numberInstances, size_t& dimension) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  dimension = 0;
  numberInstances = 0;

  if (!myfile.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }

  while (!myfile.eof()) {
    std::getline(myfile, line);
  if (line.empty())
    continue;
  if (dimension == 0) {
    dimension = std::count(line.begin(), line.end(), ',');
  } else if (dimension - std::count(line.begin(), line.end(), ',') != 0) {
    std::string msg = "Columns missing in line ";
    msg.append(std::to_string(numberInstances));
    msg.append(" in file ");
    msg.append(filename);
    throw sgpp::base::file_exception(msg.c_str());
  }
  numberInstances++;
  }
  myfile.close();
}



void CSVTools::writeNewTrainingDataEntry(const std::string& CSVLine,
    sgpp::base::DataMatrix& destination,
    size_t instanceNo) {
  size_t cur_pos = 0;
  size_t cur_find = 0;
  size_t dim = destination.getNcols();
  std::string cur_value;
  double dbl_cur_value;

  for (size_t i = 0; i < dim; i++) {
    cur_find = CSVLine.find(",", cur_pos);
    cur_value = CSVLine.substr(cur_pos, cur_find - cur_pos);
    dbl_cur_value = atof(cur_value.c_str());
    destination.set(instanceNo, i, dbl_cur_value);
    cur_pos = cur_find + 1;
  }
}

void CSVTools::writeNewClass(const std::string& CSVLine,
                              sgpp::base::DataVector& destination, size_t instanceNo) {
  size_t cur_pos = CSVLine.find_last_of(",");
  std::string cur_value = CSVLine.substr(cur_pos + 1);
  double dbl_cur_value = atof(cur_value.c_str());
  destination.set(instanceNo, dbl_cur_value);
}
}  // namespace datadriven
}  // namespace sgpp

