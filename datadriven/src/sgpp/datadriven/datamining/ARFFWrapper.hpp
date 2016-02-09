/*
 * ARFFWrapper.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef ARFFWRAPPER_HPP_
#define ARFFWRAPPER_HPP_

#include <sgpp/datadriven/datamining/DataWrapper.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>


#include <sgpp/globaldef.hpp>


#include <string>


namespace SGPP {
namespace datadriven {

class ARFFWrapper: public DataWrapper {
 public:
  ARFFWrapper(std::string filename);

  virtual ~ARFFWrapper();

  Dataset nextSamples(int how_many);

  Dataset allSamples();

protected:
  /**
   * Reads an ARFF file.
   *
   * @param filename filename of the file to be read
   * @return ARFF as Dataset
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
  static void readARFFSize(const std::string& filename, size_t& numberInstances,
                           size_t& dimension);

  static void readARFFSizeFromString(const std::string& content,
                                     size_t& numberInstances, size_t& dimension);

  /**
   * stores the attribute info of one instance into a SGPP::base::DataMatrix
   *
   * @param arffLine the string that contains the instance's values
   * @param destination SGPP::base::DataMatrix into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewTrainingDataEntry(const std::string& arffLine,
                                        SGPP::base::DataMatrix& destination, size_t instanceNo);

  /**
   * stores the class info of one instance into a SGPP::base::DataVector
   *
   * @param arffLine the string that contains the instance's class
   * @param destination SGPP::base::DataVector into which the instance is stored
   * @param instanceNo the number of the instance
   */
  static void writeNewClass(const std::string& arffLine,
                            SGPP::base::DataVector& destination, size_t instanceNo);


  size_t seed;
  size_t dimension;
  size_t numberInstances;

};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* ARFFWRAPPER_HPP_ */
