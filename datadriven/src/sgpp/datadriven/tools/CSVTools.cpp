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

#include <math.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <regex>


namespace sgpp {
namespace datadriven {

Dataset CSVTools::readCSVFromFile(const std::string& filename,
                                  bool skipFirstLine,
                                  bool hasTargets,
                                  size_t instanceCutoff,
                                  std::vector<size_t> selectedCols,
                                  std::vector<double> selectedTargets) {
  std::ifstream stream(filename.c_str());
  if (!stream.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  Dataset d =
    readCSV(stream, skipFirstLine, hasTargets, instanceCutoff, selectedCols, selectedTargets);
  stream.close();
  return d;
}

void CSVTools::readCSVSizeFromFile(const std::string& filename,
                                   size_t& numberInstances,
                                   size_t& dimension,
                                   bool skipFirstLine,
                                   bool hasTargets,
                                   std::vector<double> selectedTargets) {
  std::ifstream stream(filename.c_str());
  if (!stream.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  readCSVSize(stream, numberInstances, dimension, skipFirstLine, hasTargets, selectedTargets);
  stream.close();
}

Dataset CSVTools::readCSV(std::istream& stream,
                          bool skipFirstLine,
                          bool hasTargets,
                          size_t instanceCutoff,
                          std::vector<size_t> selectedCols,
                          std::vector<double> selectedTargets) {
  size_t maxInst = 0;
  size_t maxDim = 0;
  size_t dimension = 0;
  size_t numberInstances = 0;
  readCSVSize(stream, maxInst, maxDim, skipFirstLine, hasTargets, selectedTargets);
  // reset the stream
  stream.clear();
  stream.seekg(0, std::ios::beg);
  // make sure selectedCols has admissible values if it is not empty
  if (selectedCols.size() > 0) {
    if (*std::max_element(selectedCols.begin(), selectedCols.end()) >= maxDim) {
      throw sgpp::base::file_exception("readCSV invalid column selection");
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
  while (!stream.eof()) {
    std::getline(stream, line);
    if (skipFirstLine) {
      skipFirstLine = false;
      continue;
    }
    if (line.empty()) {
      continue;
    }
    rowEntries = tokenizeLine(line);
    // if we want a target, we remove it from line as we process it
    if (hasTargets) {
      double cl = rowEntries.back();
      rowEntries.pop_back();  // removes last element from vector
      // if no classes are specified, always accept the line
      bool isSelectedClass = selectedTargets.size() == 0;
      // TODO(sebastian): this float comp is probably implemented in SGPP
      for (size_t i = 0; i < selectedTargets.size(); i++) {
        isSelectedClass = isSelectedClass ||
          (fabs(cl - selectedTargets.at(i)) < 0.001);
      }
      if (isSelectedClass) {
        dataset.getTargets().set(rowIndex, cl);
      } else {
        // line has class which is not in selectedTargets
        continue;
      }
    }
    // We want all columns as dimensions
    if (selectedCols.size() == 0) {
      for (size_t i = 0; i < rowEntries.size(); i++) {
        dataset.getData().set(rowIndex, i, rowEntries.at(i));
      }
    } else {
      // only write dimensions specified in selectedCols
      for (size_t i = 0; i < selectedCols.size(); i++) {
        dataset.getData().set(rowIndex, i, rowEntries.at(selectedCols.at(i)));
      }
    }
    rowIndex++;
    if (rowIndex >= numberInstances) {
      break;
    }
  }
  return dataset;
}

void CSVTools::readCSVSize(std::istream& stream,
                           size_t& numberInstances,
                           size_t& dimension,
                           bool skipFirstLine,
                           bool hasTargets,
                           std::vector<double> selectedTargets) {
  std::string line;
  dimension = 0;
  numberInstances = 0;
  while (!stream.eof()) {
    std::getline(stream, line);
    if (line.empty()) {
      continue;
    }
    if (skipFirstLine) {
      skipFirstLine = false;
      continue;
    }
    // Relies on first line being a correct instance line
    if (dimension == 0) {
      dimension = std::count(line.begin(), line.end(), ',');
      dimension += hasTargets ? 0 : 1;
    } else if (dimension - std::count(line.begin(), line.end(), ',') !=
               (hasTargets ? 0 : 1)) {
      std::string msg = "readCSVSize: Columns missing in line ";
      msg.append(std::to_string(numberInstances));
      throw sgpp::base::file_exception(msg.c_str());
    }
    if (hasTargets && selectedTargets.size() > 0) {
      size_t cur_pos = line.find_last_of(",");
      std::string cur_value = line.substr(cur_pos + 1);
      double cl = atof(cur_value.c_str());
      for (size_t i = 0; i < selectedTargets.size(); i++) {
        // TODO(sebastian): Find better solution for float comp
        if (fabs(cl - selectedTargets.at(i)) < 0.001) {
          numberInstances++;
          break;
        }
      }
    } else {
      numberInstances++;
    }
  }
}

std::vector<double> CSVTools::tokenizeLine(const std::string& line) {
  sgpp::datadriven::StringTokenizer tokenizer;
  std::vector<std::string> toks;
  tokenizer.tokenize(line, ",", toks);
  std::vector<double> vals;
  for (size_t i = 0; i < toks.size(); i++) {
    vals.push_back(atof(toks.at(i).c_str()));
  }
  return vals;
}

void CSVTools::writeMatrixToCSVFile(const std::string& path, DataMatrix matrix){

 std::cout << "Writing to file "+path +".csv"<< std::endl;
 std::ofstream output;


 output.open(path+".csv");

 DataVector row(matrix.getNcols());

 for (size_t index=0; index<matrix.getNrows(); index++)
 {
   matrix.getRow(index, row);
   std::string line = row.toString();
   line.erase(line.begin());
   line.erase(line.end()-1);
   output << line+"\n";
 }

 output.close();
}

}  // namespace datadriven
}  // namespace sgpp
