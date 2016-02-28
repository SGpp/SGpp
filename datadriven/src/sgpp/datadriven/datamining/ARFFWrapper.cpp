/*
 * ARFFWrapper.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/datamining/ARFFWrapper.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <ctime>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

namespace SGPP {
namespace datadriven {

ARFFWrapper::ARFFWrapper(datadriven::DataMiningConfiguration& config)
    : DataWrapper(config), seed(0), dimension(0), numberInstances(0) {
  std::string line;
  std::ifstream myfile(filename.c_str());

  // set seed to the current time in seconds
  seed = std::time(NULL);

  if (!myfile.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw new SGPP::base::file_exception(msg.c_str());
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
    } else if (!line.empty()) {
      numberInstances++;
    }
  }

  myfile.close();
}

ARFFWrapper::~ARFFWrapper() {}

Dataset ARFFWrapper::nextSamples(int how_many) {
  // TODO(lettrich): implement
  return Dataset();
}

Dataset ARFFWrapper::allSamples() {
  Dataset dataset(numberInstances, dimension);
  readARFF(filename, dataset);
  return dataset;
}

void ARFFWrapper::readARFF(const std::string& filename, Dataset& dataset) {
  std::string line;
  std::ifstream myfile(filename.c_str());
  size_t numberInstances;
  size_t dimension;
  bool dataReached = false;
  size_t instanceNo = 0;

  readARFFSize(filename, numberInstances, dimension);
  //  Dataset dataset(numberInstances, dimension);

  while (!myfile.eof()) {
    std::getline(myfile, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (dataReached && !line.empty()) {
      writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }

    if (line.find("@DATA", 0) != line.npos) {
      dataReached = true;
    }
  }

  myfile.close();
}

void ARFFWrapper::readARFFSize(const std::string& filename, size_t& numberInstances,
                               size_t& dimension) {}

void ARFFWrapper::readARFFSizeFromString(const std::string& content, size_t& numberInstances,
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

void ARFFWrapper::readARFFFromString(const std::string& content, Dataset& dataset) {
  std::string line;
  std::stringstream contentStream;
  contentStream << content;
  size_t numberInstances;
  size_t dimension;
  bool dataReached = false;
  size_t instanceNo = 0;

  readARFFSizeFromString(content, numberInstances, dimension);
  //  Dataset dataset(numberInstances, dimension);

  while (!contentStream.eof()) {
    std::getline(contentStream, line);
    std::transform(line.begin(), line.end(), line.begin(), toupper);

    if (dataReached && !line.empty()) {
      writeNewClass(line, dataset.getTargets(), instanceNo);
      writeNewTrainingDataEntry(line, dataset.getData(), instanceNo);
      instanceNo++;
    }

    if (line.find("@DATA", 0) != line.npos) {
      dataReached = true;
    }
  }

  //  return dataset;
}

void ARFFWrapper::writeNewTrainingDataEntry(const std::string& arffLine,
                                            SGPP::base::DataMatrix& destination,
                                            size_t instanceNo) {
  size_t cur_pos = 0;
  size_t cur_find = 0;
  size_t dim = destination.getNcols();
  std::string cur_value;
  float_t dbl_cur_value;

  for (size_t i = 0; i < dim; i++) {
    cur_find = arffLine.find(",", cur_pos);
    cur_value = arffLine.substr(cur_pos, cur_find - cur_pos);
    dbl_cur_value = atof(cur_value.c_str());
    destination.set(instanceNo, i, dbl_cur_value);
    cur_pos = cur_find + 1;
  }
}

void ARFFWrapper::writeNewClass(const std::string& arffLine, SGPP::base::DataVector& destination,
                                size_t instanceNo) {
  size_t cur_pos = arffLine.find_last_of(",");
  std::string cur_value = arffLine.substr(cur_pos + 1);
  float_t dbl_cur_value = atof(cur_value.c_str());
  destination.set(instanceNo, dbl_cur_value);
}

} /* namespace datadriven */
} /* namespace SGPP */
