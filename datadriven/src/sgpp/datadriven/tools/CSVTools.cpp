// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
// Created on: 18.12.2017
// Author: Eric Koepke

#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>

#include <sgpp/globaldef.hpp>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

namespace sgpp {
namespace datadriven {

Dataset CSVTools::readCSV(const std::string& filename, bool skipFirstLine, bool hasTargets) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  size_t numberInstances;
  size_t dimension;
  size_t instanceNo = 0;

  readCSVSize(filename, numberInstances, dimension, hasTargets);
  if (skipFirstLine && !myfile.eof()) {
    numberInstances--;
    std::getline(myfile, line);
  }
  Dataset dataset(numberInstances, dimension);

  while (!myfile.eof()) {
    std::getline(myfile, line);

    if (!line.empty()) {
      if (hasTargets)
        writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }
  }

  myfile.close();

  return dataset;
}

Dataset CSVTools::readCSVPartial(const std::string& filename,
                                 bool skipFirstLine,
                                 bool hasTargets,
                                 size_t instanceCutoff,
                                 std::vector<size_t> selectedCols,
                                 std::vector<double> selectedTargets)
{
  size_t maxInst = 0;
  size_t maxDim = 0;
  size_t dimension = 0;
  size_t numberInstances = 0;
  readCSVSizePartial(filename, maxInst, maxDim, skipFirstLine,
                     hasTargets, selectedTargets);
  // make sure selectedCols has admissible values if it is not empty
  if(selectedCols.size() > 0) {
    if(*std::max_element(selectedCols.begin(), selectedCols.end()) >= maxDim) {
      throw sgpp::base::file_exception("readCSVPartial invalid col selection");
    }
    dimension = selectedCols.size();
  } else {
    dimension = maxDim;
  }
  numberInstances = std::min(maxInst, instanceCutoff);
  Dataset dataset(numberInstances, dimension);
  std::string line;
  std::ifstream myfile(filename.c_str());
  size_t rowIndex = 0;
  std::vector<double> rowEntries;
  while(!myfile.eof()) {
    std::getline(myfile, line);
    if(skipFirstLine) {
      skipFirstLine = false;
      continue;
    }
    if(!line.empty()) {
      rowEntries = tokenizeLine(line);
      // if we want a target, we remove it from line as we process it
      if(hasTargets) {
        double cl = rowEntries.back();
        rowEntries.pop_back(); // removes last element from vector
        // if no classes are specified, always accept the line
        bool isSelectedClass = selectedTargets.size() == 0;
        // TODO(sebastian): this float comp is probably implemented in SGPP
        for(size_t i = 0; i < selectedTargets.size(); i++) {
          isSelectedClass = isSelectedClass ||
            (fabs(cl - selectedTargets.at(i)) < 0.001);
        }
        if(isSelectedClass) {
          dataset.getTargets().set(rowIndex, cl);
        } else {
          // line has class which is not in selectedTargets
          continue;
        }
      }
      // We want all columns as dimensions
      if(selectedCols.size() == 0) {
        for(size_t i = 0; i < rowEntries.size(); i++) {
          dataset.getData().set(rowIndex, i, rowEntries.at(i));
        }
      } else {
        // only write dimensions specified in selectedCols
        for(size_t i = 0; i < selectedCols.size(); i++) {
          dataset.getData().set(rowIndex, i, rowEntries.at(selectedCols.at(i)));
        }
      }
      rowIndex++;
      if(rowIndex >= numberInstances) {
        break;
      }
    }
  }
  myfile.close();
  return dataset;
}

std::vector<double> CSVTools::tokenizeLine(const std::string& line) {
  sgpp::datadriven::StringTokenizer tokenizer;
  std::vector<std::string> toks;
  tokenizer.tokenize(line, ",", toks);
  std::vector<double> vals;
  for(size_t i = 0; i < toks.size(); i++) {
    vals.push_back(atof(toks.at(i).c_str()));
  }
  return vals;
}

void CSVTools::readCSVSize(const std::string& filename,
                             size_t& numberInstances, size_t& dimension, bool hasTargets) {
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
    if (!hasTargets)
      dimension++;
  } else if (dimension - std::count(line.begin(), line.end(), ',') != (hasTargets ? 0 : 1)) {
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

void CSVTools::readCSVSizePartial(const std::string& filename,
                                  size_t& numberInstances,
                                  size_t& dimension,
                                  bool skipFirstLine,
                                  bool hasTargets,
                                  std::vector<double> selectedTargets)
{
  std::string line;
  std::ifstream myfile(filename.c_str());
  dimension = 0;
  numberInstances = 0;
  if (!myfile.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  while(!myfile.eof()) {
    std::getline(myfile, line);
    if(line.empty()) {
      continue;
    }
    if(skipFirstLine) {
      skipFirstLine = false;
      continue;
    }
    if(dimension == 0) {
      dimension = std::count(line.begin(), line.end(), ',');
      dimension += hasTargets ? 0 : 1;
    } else if (dimension - std::count(line.begin(), line.end(), ',') !=
               (hasTargets ? 0 : 1)) {
      std::string msg = "Columns missing in line ";
      msg.append(std::to_string(numberInstances));
      msg.append(" in file ");
      msg.append(filename);
      throw sgpp::base::file_exception(msg.c_str());
    }
    if(hasTargets && selectedTargets.size() > 0) {
      size_t cur_pos = line.find_last_of(",");
      std::string cur_value = line.substr(cur_pos + 1);
      double cl = atof(cur_value.c_str());
      for(size_t i = 0; i < selectedTargets.size(); i++) {
        //TODO(sebastian): Find better solution for float comp
        if(fabs(cl - selectedTargets.at(i)) < 0.001) {
          numberInstances++;
          break;
        }
      }
    } else {
      numberInstances++;
    }
  }
  myfile.close();
}

void CSVTools::writeNewTrainingDataEntry(const std::string& CSVLine,
    sgpp::base::DataMatrix& destination,
    size_t instanceNo) {
  size_t dim = destination.getNcols();
  double dbl_cur_value;

  sgpp::datadriven::StringTokenizer tokenizer;
  std::vector<std::string> tokens;
  tokenizer.tokenize(CSVLine, ",", tokens);

  for (size_t i = 0; i < dim; i++) {
    dbl_cur_value = atof(tokens[i].c_str());
    destination.set(instanceNo, i, dbl_cur_value);
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

