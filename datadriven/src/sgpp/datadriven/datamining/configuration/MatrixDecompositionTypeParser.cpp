/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * MatrixDecompositionTypeParser.cpp
 *
 * Created on: Jan 30, 2018
 *     Author: Kilian Röhner
 */

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

MatrixDecompositionType MatrixDecompositionTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower.compare("cg") == 0) {
    return sgpp::datadriven::MatrixDecompositionType::LU;
  } else if (inputLower.compare("eigen") == 0) {
    return sgpp::datadriven::MatrixDecompositionType::Eigen;
  } else if (inputLower.compare("chol") == 0) {
    return sgpp::datadriven::MatrixDecompositionType::Chol;
  } else if (inputLower.compare("denseichol") == 0) {
    return sgpp::datadriven::MatrixDecompositionType::DenseIchol;
  } else if (inputLower.compare("orthoadapt") == 0) {
    return sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any "
        "known MatrixDecompositionType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string &MatrixDecompositionTypeParser::toString(MatrixDecompositionType type) {
  return matrixDecompositionTypeMap.at(type);
}

const MatrixDecompositionTypeParser::MatrixDecompositionTypeMap_t
    MatrixDecompositionTypeParser::matrixDecompositionTypeMap = []() {
  return MatrixDecompositionTypeMap_t{
      std::make_pair(MatrixDecompositionType::LU, "LU"),
      std::make_pair(MatrixDecompositionType::Eigen, "Eigen"),
      std::make_pair(MatrixDecompositionType::Chol, "Chol"),
      std::make_pair(MatrixDecompositionType::DenseIchol, "DenseIchol"),
      std::make_pair(MatrixDecompositionType::OrthoAdapt, "OrthoAdapt")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
