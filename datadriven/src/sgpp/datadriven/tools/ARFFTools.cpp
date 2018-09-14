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

Dataset ARFFTools::readARFF(const std::string& filename, bool hasTargets) {
  // TODO(fuchsgruber): No idea if this arff interface really can handle data without classes
  std::string line;
  std::ifstream myfile(filename.c_str());
  if (!myfile) {
    const auto msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  size_t numberInstances;
  size_t dimension;
  bool dataReached = false;
  size_t instanceNo = 0;

  readARFFSize(filename, numberInstances, dimension);
  Dataset dataset(numberInstances, dimension);

  while (!myfile.eof()) {
    std::getline(myfile, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (dataReached && !line.empty()) {
      if (hasTargets)
        writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }

    if (line.find("@DATA", 0) != line.npos) {
      dataReached = true;
    }
  }

  myfile.close();

  return dataset;
}


Dataset ARFFTools::readARFFPartial(const std::string& filename,
                                   bool hasTargets,
                                   size_t instanceCutoff,
                                   std::vector<size_t> selectedCols,
                                   std::vector<double> selectedTargets) {
  size_t maxInst = 0;
  size_t maxDim = 0;
  size_t dimension = 0;
  size_t numberInstances = 0;
  readARFFSizePartial(filename, maxInst, maxDim, hasTargets, selectedTargets);
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
  myfile.close();
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

void ARFFTools::readARFFSize(const std::string& filename, size_t& numberInstances,
                             size_t& dimension) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  dimension = 0;
  numberInstances = 0;

  if (!myfile) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }

  while (!myfile.eof()) {
    std::getline(myfile, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (line.find("@ATTRIBUTE class", 0) != line.npos) {
    } else if (line.find("@ATTRIBUTE CLASS", 0) != line.npos) {
    } else if (line.find("@ATTRIBUTE", 0) != line.npos) {
      dimension++;
    } else if (line.find("@DATA", 0) != line.npos) {
      numberInstances = 0;
    } else if (line.find("% DATA SET SIZE ", 0) != line.npos) {
      numberInstances = std::stoi(line.substr(strlen("% DATA SET SIZE ")));
      std::cout << "Set number instances from comment to " << numberInstances << std::endl;
      break;
    } else if (!line.empty()) {
      numberInstances++;
    }
  }

  myfile.close();
}

void ARFFTools::readARFFSizePartial(const std::string& filename,
                                           size_t& numberInstances,
                                           size_t& dimension,
                                           bool hasTargets,
                                           std::vector<double> selectedTargets) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  dimension = 0;
  numberInstances = 0;
  if (!myfile) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  while (!myfile.eof()) {
    std::getline(myfile, line);
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
void ARFFTools::readARFFSizeFromString(const std::string& content, size_t& numberInstances,
                                       size_t& dimension) {
  std::string line;
  std::istringstream contentStream(content);

  dimension = 0;
  numberInstances = 0;

  while (!contentStream.eof()) {
    std::getline(contentStream, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (line.find("@ATTRIBUTE class", 0) != line.npos) {
    } else if (line.find("@ATTRIBUTE CLASS", 0) != line.npos) {
    } else if (line.find("@ATTRIBUTE", 0) != line.npos) {
      dimension++;
    } else if (line.find("@DATA", 0) != line.npos) {
      numberInstances = 0;
    } else if (!line.empty()) {
      numberInstances++;
    }
  }
}

Dataset ARFFTools::readARFFFromString(const std::string& content, bool hasTargets) {
  std::string line;
  std::stringstream contentStream;
  contentStream << content;
  size_t numberInstances;
  size_t dimension;
  bool dataReached = false;
  size_t instanceNo = 0;

  ARFFTools::readARFFSizeFromString(content, numberInstances, dimension);
  Dataset dataset(numberInstances, dimension);

  while (!contentStream.eof()) {
    std::getline(contentStream, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (dataReached && !line.empty()) {
      if (hasTargets)
        writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }

    if (line.find("@DATA", 0) != line.npos) {
      dataReached = true;
    }
  }

  return dataset;
}

void ARFFTools::writeNewTrainingDataEntry(const std::string& arffLine,
                                          sgpp::base::DataMatrix& destination, size_t instanceNo) {
  size_t cur_pos = 0;
  size_t cur_find = 0;
  size_t dim = destination.getNcols();
  std::string cur_value;
  double dbl_cur_value;

  for (size_t i = 0; i < dim; i++) {
    cur_find = arffLine.find(",", cur_pos);
    cur_value = arffLine.substr(cur_pos, cur_find - cur_pos);
    dbl_cur_value = atof(cur_value.c_str());
    destination.set(instanceNo, i, dbl_cur_value);
    cur_pos = cur_find + 1;
  }
}

void ARFFTools::writeNewClass(const std::string& arffLine, sgpp::base::DataVector& destination,
                              size_t instanceNo) {
  size_t cur_pos = arffLine.find_last_of(",");
  std::string cur_value = arffLine.substr(cur_pos + 1);
  double dbl_cur_value = atof(cur_value.c_str());
  destination.set(instanceNo, dbl_cur_value);
}

// void ARFFTools::writeAlpha(std::string tfilename, sgpp::base::DataVector& source)
// {
//
// }

// void ARFFTools::readAlpha(std::string tfilename, sgpp::base::DataVector& destination)
// {
//
// }

}  // namespace datadriven
}  // namespace sgpp
