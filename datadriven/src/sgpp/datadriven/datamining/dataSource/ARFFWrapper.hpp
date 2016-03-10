/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ARFFWrapper.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef ARFFWRAPPER_HPP_
#define ARFFWRAPPER_HPP_

#include <string>

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/datadriven/datamining/dataSource/DataWrapper.hpp>

namespace sgpp {
namespace datadriven {

class ARFFWrapper : public DataWrapper {
 public:
  ARFFWrapper(datadriven::DataMiningConfigJsonParser& config);

  virtual ~ARFFWrapper();

  Dataset nextSamples(int how_many);

  Dataset allSamples();

 protected:
  /**
   * Reads an ARFF file.
   *
   * @param filename filename of the file to be read
   * @param[out] dataset ARFF as Dataset
   */
  static void readARFF(const std::string& filename, Dataset& dataset);

  static void readARFFFromString(const std::string& content, Dataset& dataset);

  /**
   * Reads the size of an ARFF file.
   *
   * @param filename filename of the file to be read
   * @param[out] numberInstances number of instances in the dataset
   * @param[out] dimension number of dimensions in the dataset
   */
  static void readARFFSize(const std::string& filename, size_t& numberInstances, size_t& dimension);

  static void readARFFSizeFromString(const std::string& content, size_t& numberInstances,
                                     size_t& dimension);

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
  static void writeNewClass(const std::string& arffLine, sgpp::base::DataVector& destination,
                            size_t instanceNo);

  size_t seed;
  size_t dimension;
  size_t numberInstances;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* ARFFWRAPPER_HPP_ */
