// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class MatrixDecompositionTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::datadriven::MatrixDecompositionType.
   * Throws if there is no valid representation
   * @param input case insensitive string representation of a
   * #sgpp::datadriven::MatrixDecompositionType.
   * @return the corresponding #sgpp::datadriven::MatrixDecompositionType.
   */
  static MatrixDecompositionType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::datadriven::MatrixDecompositionType.
   * @param type enum value.
   * @return string representation of a #sgpp::datadriven::MatrixDecompositionType.
   */
  static const std::string &toString(MatrixDecompositionType type);

 private:
  typedef std::map<MatrixDecompositionType, std::string> MatrixDecompositionTypeMap_t;

  /**
   * Map containing all values of  #sgpp::datadriven::MatrixDecompositionType and
   * the corresponding string representation.
   */
  static const MatrixDecompositionTypeMap_t matrixDecompositionTypeMap;
};
} /* namespace datadriven */
} /* namespace sgpp */
