// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/CSVWriter.hpp>
#include <sgpp/base/exception/file_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <math.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

void CSVWriter::writeCSVToFile(const std::string& filename, const Dataset& dataset, bool hasTargets,
                               size_t instanceCutoff, std::vector<size_t> selectedCols) {
  std::ofstream stream(filename.c_str());
  if (!stream.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    throw sgpp::base::file_exception(msg.c_str());
  }
  writeCSV(stream, dataset, hasTargets, instanceCutoff, selectedCols);
  stream.close();
}

void CSVWriter::writeCSV(std::ostream& stream, const Dataset& dataset, bool hasTargets,
                         size_t instanceCutoff, std::vector<size_t> selectedCols) {
  size_t maxInst = dataset.getNumberInstances();
  size_t maxDim = dataset.getDimension();
  size_t dimension = 0;
  size_t numberInstances = 0;

  // reset the stream
  stream.clear();
  // make sure selectedCols has admissible values if it is not empty
  if (selectedCols.size() > 0) {
    if (*std::max_element(selectedCols.begin(), selectedCols.end()) >= maxDim) {
      throw sgpp::base::file_exception("writeCSV invalid column selection");
    }
    dimension = selectedCols.size();
  } else {
    dimension = maxDim;
  }
  numberInstances = std::min(maxInst, instanceCutoff);

  auto& m = dataset.getData();
  auto& t = dataset.getTargets();

  std::string line;
  size_t rowIndex = 0;

  for (size_t i = 0; i < m.getNrows(); ++i) {
    sgpp::base::DataVector rowEntries(dimension);

    // We want all columns as dimensions
    if (selectedCols.size() == 0) {
      for (size_t j = 0; j < maxDim; ++j) {
        rowEntries.set(j, m.get(i, j));
      }
    } else {
      // only write dimensions specified in selectedCols
      for (size_t j = 0; j < selectedCols.size(); ++j) {
        rowEntries.set(j, m.get(i, selectedCols[j]));
      }
    }

    line = stringifyLine(rowEntries);
    stream << line;
    if (hasTargets) stream << "," << std::to_string(t.get(i));
    stream << std::endl;

    rowIndex++;
    if (rowIndex >= numberInstances) {
      break;
    }
  }
}

std::string CSVWriter::stringifyLine(const sgpp::base::DataVector& values) {
  std::string line;
  for (double val : values) {
    line.append(std::to_string(val));
    line.append(",");
  }
  line.pop_back();
  return line;
}

}  // namespace datadriven
}  // namespace sgpp
