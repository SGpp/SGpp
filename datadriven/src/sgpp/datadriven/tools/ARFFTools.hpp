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
   * Reads the size of a ARFF file.
   *
   * @param stream contains the raw data
   * @param[out] numberInstances number of instances in the dataset
   * @param[out] number of columns (dimensions) in the dataset
   * @param hasTargets whether the csv has a columns for targets
   *        (supervised learning). If true, dimension = #columns - 1
   * @param selectedTargets (see readCSVPartial). If this vector is not empty (default)
   *        numberIstances reflects only the number of instance which are admissible
   *        with respect to selectedTargets.
   *        If empty all targets are admissible and all rows are considered as instances.
   */
   static void readARFFSize(std::istream& stream,
                            size_t& numberInstances,
                            size_t& dimension,
                            bool hasTargets,
                            std::vector<double> selectedTargets);

  /**
   * Sequentially reads a ARFF file.
   *
   * @param stream contains the raw data
   * @param hasTargets whether the csv has columns for targets
   *        (supervised learning)
   * @param instanceCutoff maximal number of instances to include in the
   *        returned Dataset. May not be reached if there are less than instanceCutoff
   *        (valid w.r.t. selectedTargets) rows in csv file. If the value is -1
   *        i.e. the maximal value of the (unsigned) size_t, all valid rows are included
   * @param selectedCols which columns are written to the DataMatrix as
   *        dimensions. Order matters, i.e. selectedCols = [0, 3, 2]
   *        will result in a DataMatrix with dim0 = row0, dim1 = row3,
   *        dim2 = row2. If hasTargets=true, the last row must not be
   *        specified here as it is written to the target vector, not the
   *        DataMatrix. If empty (default) all rows (except possible the
   *        target row) will be used in ascending order.
   * @param selectedTargets filter for targets. Only applicable if
   *        hasTargets=true. Only rows with target-entry (last column)
   *        equal to one of the entries in selectedTarget are written to
   *        the dataset. All other rows are skipped. Float-comparison uses
   *        0.001 precision. This parameter is intended to use for
   *        classification with integer values as classes.
   *        If empty (default) all targets are admissible and all rows are
   *        written as rows (instances) in the dataset.
   * @return ARFF as Dataset
   */
   static Dataset readARFF(std::istream& stream,
                           bool hasTargets = true,
                           size_t instanceCutoff = -1,
                           std::vector<size_t> selectedCols = std::vector<size_t>(),
                           std::vector<double> selectedTargets = std::vector<double>());

   /**
    * Wrapper from input type: File. See readARFFSize for more details
    */
   static void readARFFSizeFromFile(const std::string& filename,
                                    size_t& numberInstances,
                                    size_t& dimension,
                                    bool hasTargets = true,
                                    std::vector<double> selectedTargets = std::vector<double>());

   /**
    * Wrapper from input type: String. See readARFFSize for more details
    */
   static void readARFFSizeFromString(const std::string& content,
                                      size_t& numberInstances,
                                      size_t& dimension,
                                      bool hasTargets = true,
                                      std::vector<double> selectedTargets = std::vector<double>());

   /**
    * Wrapper from input type: File. See readARFF for more details
    */
   static Dataset readARFFFromFile(const std::string& filename,
                                   bool hasTargets = true,
                                   size_t instanceCutoff = -1,
                                   std::vector<size_t> selectedCols = std::vector<size_t>(),
                                   std::vector<double> selectedTargets = std::vector<double>());

   /**
    * Wrapper from input type: String. See readARFF for more details
    */
   static Dataset readARFFFromString(const std::string& content,
                                     bool hasTargets = true,
                                     size_t instanceCutoff = -1,
                                     std::vector<size_t> selectedCols = std::vector<size_t>(),
                                     std::vector<double> selectedTargets = std::vector<double>());


 private:
   /**
   * Take a comma-sererated line and return its values as double values
   */
  static std::vector<double> tokenizeLine(const std::string& line);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ARFFTOOLS_HPP */
