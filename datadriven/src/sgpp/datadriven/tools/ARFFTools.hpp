// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ARFFTOOLS_HPP
#define ARFFTOOLS_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Class that provides functionality to read ARFF files.
 */
class ARFFTools {
 public:
  /**
   * Reads an ARFF file.
   *
   * @param filename filename of the file to be read
   * @param hasTargets whether the file has targest (i.e. supervised learning)
   * @return ARFF as Dataset
   */
  static Dataset readARFF(const std::string& filename, bool hasTargets = true);

  static Dataset readARFFPartial(const std::string& filename,
                                 bool hasTargets = true,
                                 size_t instanceCutoff = -1,
                                 std::vector<size_t> selectedCols = std::vector<size_t>(),
                                 std::vector<double> selectedTargets = std::vector<double>());

  /**
   * Reads an ARFF file content.
   *
   * @param content the file content to read
   * @param hasTargets whether the file has targest (i.e. supervised learning)
   * @return ARFF as Dataset
   */
  static Dataset readARFFFromString(const std::string& content, bool hasTargets = true);

  /**
   * Reads the size of an ARFF file.
   *
   * @param filename filename of the file to be read
   * @param[out] numberInstances number of instances in the dataset
   * @param[out] dimension number of dimensions in the dataset
   */
  static void readARFFSize(const std::string& filename, size_t& numberInstances,
                           size_t& dimension);

  static void readARFFSizePartial(const std::string& filename,
                                  size_t& numberInstances,
                                  size_t& dimension,
                                  bool hasTargets = true,
                                  std::vector<double> selectedTargets = std::vector<double>());

  static void readARFFSizeFromString(const std::string& content,
                                     size_t& numberInstances, size_t& dimension);



 private:
  /**
   * stores the attribute info of one instance into a sgpp::base::DataMatrix
   *
   * @param arffLine the string that contains the instance's values
   * @param destination sgpp::base::DataMatrix into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewTrainingDataEntry(const std::string& arffLine,
                                        sgpp::base::DataMatrix& destination, size_t instanceNo);

  /**
   * stores the class info of one instance into a sgpp::base::DataVector
   *
   * @param arffLine the string that contains the instance's class
   * @param destination sgpp::base::DataVector into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewClass(const std::string& arffLine,
                            sgpp::base::DataVector& destination, size_t instanceNo);

  static std::vector<double> tokenizeLine(const std::string& line);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ARFFTOOLS_HPP */
