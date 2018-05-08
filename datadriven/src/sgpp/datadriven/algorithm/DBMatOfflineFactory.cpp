/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineFactory.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::GridType;
using sgpp::base::factory_exception;

DBMatOffline* DBMatOfflineFactory::buildOfflineObject(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::base::AdpativityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  auto type = densityEstimationConfig.decomposition_;

  switch (type) {
    case (MatrixDecompositionType::Eigen):
#ifdef USE_GSL
      return new DBMatOfflineEigen();
#else
      throw factory_exception("built without GSL");
#endif /* USE_GSL */
      break;

    case (MatrixDecompositionType::LU):
#ifdef USE_GSL
      return new DBMatOfflineLU();
#else
      throw factory_exception("built without GSL");
#endif /* USE_GSL */
      break;

    case (MatrixDecompositionType::Chol):
#ifdef USE_GSL
      return new DBMatOfflineChol();
#else
      throw factory_exception("built without GSL");
#endif /* USE_GSL */
      break;

    case (MatrixDecompositionType::DenseIchol):
      return new DBMatOfflineDenseIChol();
      break;

    case (MatrixDecompositionType::OrthoAdapt):
#ifdef USE_GSL
      return new DBMatOfflineOrthoAdapt();
      break;
#else
      throw factory_exception("built without GSL");
#endif /* USE_GSL */

    default:
      throw factory_exception("Trying to build offline object from unknown decomposition type");
      return nullptr;
  }
}

DBMatOffline* DBMatOfflineFactory::buildFromFile(const std::string& fileName) {
#ifdef USE_GSL
  std::ifstream file(fileName, std::istream::in);

  if (!file) {
    throw factory_exception("Failed to open File");
  }

  std::string str;
  std::getline(file, str);
  file.close();

  std::cout << str;
  std::vector<std::string> tokens;
  StringTokenizer::tokenize(str, ",", tokens);
  std::cout << "tokens: ";
  for (auto& item : tokens) {
    std::cout << item << ",";
  }
  std::cout << std::endl;

  MatrixDecompositionType type = static_cast<MatrixDecompositionType>(std::stoi(tokens[5]));

  std::cout << "type: " << static_cast<int>(type) << std::endl;

  switch (type) {
    case (MatrixDecompositionType::Eigen):
      return new DBMatOfflineEigen(fileName);
      break;
    case (MatrixDecompositionType::LU):
      return new DBMatOfflineLU(fileName);
      break;
    case (MatrixDecompositionType::Chol):
      return new DBMatOfflineChol(fileName);
      break;
    case (MatrixDecompositionType::DenseIchol):
      return new DBMatOfflineDenseIChol(fileName);
      break;
    case (MatrixDecompositionType::OrthoAdapt):
      return new DBMatOfflineOrthoAdapt(fileName);
      break;
    default:
      throw factory_exception("Trying to build offline object from unknown decomposition type");
      return nullptr;
  }
#else
  throw factory_exception("built withot GSL");
#endif /* USE_GSL */
}

} /* namespace datadriven */
} /* namespace sgpp */
