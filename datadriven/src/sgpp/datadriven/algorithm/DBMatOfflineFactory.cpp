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

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineIChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/datamining/base/StringTokenizer.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::GridType;
using sgpp::base::algorithm_exception;

DBMatOffline* DBMatOfflineFactory::buildOfflineObject(const DBMatDensityConfiguration& config) {
  auto type = config.decomp_type_;

  switch (type) {
    case (DBMatDecompostionType::DBMatDecompEigen):
      return new DBMatOfflineEigen(config);
      break;
    case (DBMatDecompostionType::DBMatDecompLU):
      return new DBMatOfflineLU(config);
      break;
    case (DBMatDecompostionType::DBMatDecompChol):
      return new DBMatOfflineChol(config);
      break;
    case (DBMatDecompostionType::DBMatDecompIChol):
      return new DBMatOfflineIChol(config);
      break;
    default:
      throw algorithm_exception("Trying to build offline object from unknown decomposition type");
      return nullptr;
  }
}

DBMatOffline* DBMatOfflineFactory::buildFromFile(const std::string& fileName) {
  std::ifstream file(fileName, std::istream::in);

  if (!file) {
    throw algorithm_exception("Failed to open File");
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
    case (DBMatDecompostionType::DBMatDecompEigen):
      return new DBMatOfflineEigen(fileName);
      break;
    case (DBMatDecompostionType::DBMatDecompLU):
      return new DBMatOfflineLU(fileName);
      break;
    case (DBMatDecompostionType::DBMatDecompChol):
      return new DBMatOfflineChol(fileName);
      break;
    case (DBMatDecompostionType::DBMatDecompIChol):
      return new DBMatOfflineIChol(fileName);
      break;
    default:
      throw algorithm_exception("Trying to build offline object from unknown decomposition type");
      return nullptr;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
