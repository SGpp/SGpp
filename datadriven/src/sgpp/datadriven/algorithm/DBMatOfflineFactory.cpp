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
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::GridType;
using sgpp::base::factory_exception;

DBMatOffline* DBMatOfflineFactory::buildOfflineObject(const DBMatDensityConfiguration& config) {
  auto type = config.decomp_type_;

  switch (type) {
    case (DBMatDecompostionType::Eigen):
#ifdef USE_GSL
      return new DBMatOfflineEigen(config);
#else
      throw factory_exception("built withot GSL");
#endif /* USE_GSL */
      break;
    case (DBMatDecompostionType::LU):
#ifdef USE_GSL
      return new DBMatOfflineLU(config);
#else
      throw factory_exception("built withot GSL");
#endif /* USE_GSL */
      break;
    case (DBMatDecompostionType::Chol):
#ifdef USE_GSL
      return new DBMatOfflineChol(config);
#else
      throw factory_exception("built withot GSL");
#endif /* USE_GSL */
      break;
    case (DBMatDecompostionType::DenseIchol):
      return new DBMatOfflineDenseIChol(config);
      break;
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

  DBMatDecompostionType type = static_cast<DBMatDecompostionType>(std::stoi(tokens.back()));

  std::cout << "type: " << static_cast<int>(type) << std::endl;

  switch (type) {
    case (DBMatDecompostionType::Eigen):
      return new DBMatOfflineEigen(fileName);
      break;
    case (DBMatDecompostionType::LU):
      return new DBMatOfflineLU(fileName);
      break;
    case (DBMatDecompostionType::Chol):
      return new DBMatOfflineChol(fileName);
      break;
    case (DBMatDecompostionType::DenseIchol):
      return new DBMatOfflineDenseIChol(fileName);
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
