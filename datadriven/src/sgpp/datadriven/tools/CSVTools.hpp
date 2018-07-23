// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
// Created on: 18.12.2017
// Author: Eric Koepke

#ifndef CSVTOOLS_HPP
#define CSVTOOLS_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Class that provides functionality to read CSV files.
 */
class CSVTools {
 public:
  /**
   * Reads an CSV file.
   *
   * @param filename filename of the file to be read
   * @param skipFirstLine whether to skip the first line while parsing
   * @param hasTargets whether the csv has columns for targets (supervised learning)
   * @return CSV as Dataset
   */
  static Dataset readCSV(const std::string& filename, bool skipFirstLine = false,
      bool hasTargets = true);

  /**
   * Reads the size of an CSV file.
   *
   * @param filename filename of the file to be read
   * @param[out] numberInstances number of instances in the dataset
   * @param[out] dimension number of dimensions in the dataset
   * @param hasTargets whether the csv has a columns for targets (supervised learning)
   */
  static void readCSVSize(const std::string& filename, size_t& numberInstances,
                           size_t& dimension, bool hasTargets = true);

 private:
  /**
   * stores the attribute info of one instance into a sgpp::base::DataMatrix
   *
   * @param CSVLine the string that contains the instance's values
   * @param destination sgpp::base::DataMatrix into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewTrainingDataEntry(const std::string& CSVLine,
                                        sgpp::base::DataMatrix& destination, size_t instanceNo);

  /**
   * stores the class info of one instance into a sgpp::base::DataVector
   *
   * @param CSVLine the string that contains the instance's class
   * @param destination sgpp::base::DataVector into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewClass(const std::string& CSVLine,
                            sgpp::base::DataVector& destination, size_t instanceNo);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* CSVTOOLS_HPP */
