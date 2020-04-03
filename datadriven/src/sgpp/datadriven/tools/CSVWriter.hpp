// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Class that provides functionality to write a Dataset to a CSV file.
 */
class CSVWriter {
 public:
  static void writeCSV(std::ostream& stream, const Dataset& dataset, bool hasTargets = true,
                       size_t instanceCutoff = -1,
                       std::vector<size_t> selectedCols = std::vector<size_t>());

  static void writeCSVToFile(const std::string& filename, const Dataset& dataset,
                             bool hasTargets = true, size_t instanceCutoff = -1,
                             std::vector<size_t> selectedCols = std::vector<size_t>());

 private:
  /**
   * Take a vector of doubles and return a comma-separated line
   */
  static std::string stringifyLine(const sgpp::base::DataVector& values);
};

}  // namespace datadriven
}  // namespace sgpp
