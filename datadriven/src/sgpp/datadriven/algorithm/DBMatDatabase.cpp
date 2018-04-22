// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * DBMatDatabase.cpp
 *
 *  Created on: Apr 22, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

DBMatDatabase::DBMatDatabase(const std::string& filepath) {
  databaseRoot = std::make_unique<json::JSON>(filepath);
  // Get the root node of the database (list)
  if (databaseRoot->contains("database")) {
    database = (json::ListNode*)(&(*databaseRoot)["database"]);
  } else {
    std::cout << "DBMatDatabase: json database is ill formated (does not contain key \"database\")!"
        << std::endl;
  }
}

DBMatOffline* DBMatDatabase::dataMatrixFromDatabase(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdpativityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {

  try {
    std::string matrixfile = matrixFileByConfiguration(gridConfig, adaptivityConfig,
        regularizationConfig, densityEstimationConfig);
    // Offline object from file
    try {
      return sgpp::datadriven::DBMatOfflineFactory::buildFromFile(matrixfile);
    } catch (sgpp::base::factory_exception e) {
      std::cout << "DBMatOfflineFactory threw an exception while deserializing matrix "
          "decomposition: " << e.what() << std::endl;
      return nullptr;
    }

  } catch (sgpp::base::factory_exception& e) {
    // No decomposition in database
    return nullptr;
  }
}


bool DBMatDatabase::gridConfigurationMatches(json::DictNode *node,
      sgpp::base::GeneralGridConfiguration& gridConfig, size_t entry_num) {
  // Check if grid general type matches
  if (node->contains("type")) {
    // Parse the grid general type
    sgpp::base::GeneralGridType gridType;
    std::string strGridType = (*node)["type"].get();
    std::transform(strGridType.begin(), strGridType.end(), strGridType.begin(), ::tolower);
    if (strGridType == "regular") {
      gridType = sgpp::base::GeneralGridType::RegularSparseGrid;
    } else if (strGridType == "refinedcoarsened") {
      gridType = sgpp::base::GeneralGridType::RefinedCoarsenedSparseGrid;
    } else if (strGridType == "withinteractions") {
      gridType = sgpp::base::GeneralGridType::SparseGridWithInteractions;
    } else if (strGridType == "combi") {
      gridType = sgpp::base::GeneralGridType::CombiGrid;
    } else {
      std::cout << "DBMatDatabase: database entry # " << entry_num <<
              ": \"gridConfiguration\" contains invalid contain \"type\" value \"" << strGridType
              << "and therefore is ignored!" << std::endl;
      return false;
    }
    // Check if the general grid type matches
    if (gridConfig.generalType_ != gridType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"gridConfiguration\" node does not contain \"type\" key and therefore is ignored!"
        << std::endl;
    return false;
  }
  // Check if grid dimensionality matches
  if (node->contains("dimension")) {
    uint64_t gridDimension = (*node)["dimension"].getUInt();
    if (gridConfig.dim_ != gridDimension) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"gridConfiguration\" node does not contain \"dimension\" key and therefore is "
        << "ignored!" << std::endl;
    return false;
  }
  // Check if the grid level matches
  if (node->contains("level")) {
    int64_t gridLevel = (*node)["level"].getInt();
    if (gridConfig.level_ != gridLevel) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"gridConfiguration\" node does not contain \"level\" key and therefore is ignored!"
        << std::endl;
    return false;
  }
  // All relevant attributes of the general grid configuration match
  return true;
}

bool DBMatDatabase::regularizationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig, size_t entry_num) {
  // Check if the reguluarization strength matches
  if (node->contains("lambda")) {
    double lambda = (*node)["lambda"].getDouble();
    if (regularizationConfig.lambda_ != lambda) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"regularizationConfiguration\" node does not contain \"lambda\" key and therefore is "
        << "ignored!" << std::endl;
    return false;
  }
  // Regularization configuration matches
  return true;
}

bool DBMatDatabase::densityEstimationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig, size_t entry_num) {
  // Check if the decomposition type matches
  if (node->contains("decomposition")) {
    // Parse decomposition type
    sgpp::datadriven::MatrixDecompositionTypeParser typeParser;
    std::string strDecompType = (*node)["decomposition"].get();
    sgpp::datadriven::MatrixDecompositionType decompositionType = typeParser.parse(strDecompType);
    if (densityEstimationConfig.decomposition_ != decompositionType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"densityEstimationConfiguration\" node does not contain \"decomposition\" key and " <<
        "therefore is ignored!" << std::endl;
    return false;
  }
  // Density estimation configuration matches
  return true;
}

std::string DBMatDatabase::matrixFileByConfiguration(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdpativityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // Scan the entire database
  for (size_t i = 0; i < database->size(); i++) {
    json::DictNode* entry = (json::DictNode*)(&((*database)[i]));
    // Check if the entry matches the grid configuration
    if (entry->contains("gridConfiguration")) {
      json::DictNode *gridConfigNode = (json::DictNode*)(&(*entry)["gridConfiguration"]);
      if (!gridConfigurationMatches(gridConfigNode, gridConfig, i)) continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"gridConfiguration\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // Check if the entry matches the regularization configuration
    if (entry->contains("regularizationConfiguration")) {
      json::DictNode *regularizationConfigNode =
          (json::DictNode*)(&(*entry)["regularizationConfiguration"]);
      if (!regularizationConfigurationMatches(regularizationConfigNode, regularizationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"regularizationConfiguration\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // Check if the entry matches the density estimation configuration
    if (entry->contains("densityEstimationConfiguration")) {
      json::DictNode *densityEstimationConfiNode =
          (json::DictNode*)(&(*entry)["densityEstimationConfiguration"]);
      if (!densityEstimationConfigurationMatches(densityEstimationConfiNode,
          densityEstimationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"densityEstimationConfiguration\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // All three configurations match -> return this entry
    if (entry->contains("filepath")) {
      return ((*entry)["filepath"]).get();
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " matches but does not contain a " <<
          "\"filepath\" key and therefore is ignored!" << std::endl;
    }
  }
  throw sgpp::base::factory_exception("No decomposition for given configuration in database");
}
} /* namespace datadriven */
} /* namespace sgpp */




