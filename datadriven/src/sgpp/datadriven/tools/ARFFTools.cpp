// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>

namespace sgpp {
namespace datadriven {


void ARFFTools::readARFFSizeFromFile(const std::string& filename,
                             size_t& numberInstances,
                             size_t& dimension,
                             bool hasTargets,
                             std::vector<double> selectedTargets)
{
  std::ifstream stream(filename.c_str());
  if (!stream) {
    std::string msg = "readARFFSizeFromFile: Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  readARFFSize(stream, numberInstances, dimension, hasTargets, selectedTargets);
  stream.close();
}

void ARFFTools::readARFFSizeFromString(const std::string& content,
                                       size_t& numberInstances,
                                       size_t& dimension,
                                       bool hasTargets,
                                       std::vector<double> selectedTargets) {
  std::istringstream stream(content);
  readARFFSize(stream, numberInstances, dimension, hasTargets, selectedTargets);
}

Dataset ARFFTools::readARFFFromFile(const std::string& filename,
                                    bool hasTargets,
                                    size_t instanceCutoff,
                                    std::vector<size_t> selectedCols,
                                    std::vector<double> selectedTargets) {
  std::ifstream stream(filename.c_str());
  if (!stream) {
    std::string msg = "readARFFFromFile: Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  Dataset d = readARFF(stream, hasTargets, instanceCutoff, selectedCols, selectedTargets);
  stream.close();
  return d;
}

Dataset ARFFTools::readARFFFromString(const std::string& content,
                                      bool hasTargets,
                                      size_t instanceCutoff,
                                      std::vector<size_t> selectedCols,
                                      std::vector<double> selectedTargets) {
  std::istringstream stream(content);
  return readARFF(stream, hasTargets, instanceCutoff, selectedCols, selectedTargets);
}

void ARFFTools::readARFFSize(std::istream& stream,
                             size_t& numberInstances,
                             size_t& dimension,
                             bool hasTargets,
                             std::vector<double> selectedTargets)
{
  std::string line;
  dimension = 0;
  numberInstances = 0;
  while (!stream.eof()) {
    std::getline(stream, line);
    // We don't care about the attribute specification. Just skip it
    if (line.find("%", 0) != line.npos || line.find("@", 0) != line.npos) {
      continue;
    }
    if(line.empty()) {
      continue;
    }
    // Relies on first line being a correct instance line
    if(dimension == 0) {
      dimension = std::count(line.begin(), line.end(), ',');
      dimension += hasTargets ? 0 : 1;
    } else if (dimension - std::count(line.begin(), line.end(), ',') !=
               (hasTargets ? 0 : 1)) {
      std::string msg = "readARFFSize: Columns missing in line ";
      msg.append(std::to_string(numberInstances));
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
}

Dataset ARFFTools::readARFF(std::istream& stream,
                            bool hasTargets,
                            size_t instanceCutoff,
                            std::vector<size_t> selectedCols,
                            std::vector<double> selectedTargets) {
  size_t maxInst = 0;
  size_t maxDim = 0;
  size_t dimension = 0;
  size_t numberInstances = 0;
  readARFFSize(stream, maxInst, maxDim, hasTargets, selectedTargets);
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
  size_t rowIndex = 0;
  std::vector<double> rowEntries;
  while(!stream.eof()) {
    std::getline(stream, line);
    if (line.find("%", 0) != line.npos || line.find("@", 0) != line.npos) {
      continue;
    }
    if(line.empty()) {
      continue;
    }
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
  return dataset;
}

std::vector<double> ARFFTools::tokenizeLine(const std::string& line) {
  sgpp::datadriven::StringTokenizer tokenizer;
  std::vector<std::string> toks;
  tokenizer.tokenize(line, ",", toks);
  std::vector<double> vals;
  for(size_t i = 0; i < toks.size(); i++) {
    vals.push_back(atof(toks.at(i).c_str()));
  }
  return vals;
}

}  // namespace datadriven
}  // namespace sgpp
