// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GeneralGridTypeParser.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

// Globally define all keys in a database entry
const std::string keyGridConfiguration = "gridConfiguration";
const std::string keyRegularizationConfiguration = "regularizationConfiguration";
const std::string keyDensityEstimationConfiguration = "densityEstimationConfiguration";
const std::string keyGridType = "type";
const std::string keyGridDimension = "dimension";
const std::string keyGridLevel = "level";
const std::string keyRegularizationStrength = "lambda";
const std::string keyDecompositionType = "decomposition";
const std::string keyFilepath = "filepath";

DBMatDatabase::DBMatDatabase(const std::string& filepath) {
  databaseFilepath = filepath;
  databaseRoot = std::make_unique<json::JSON>(filepath);
  // Get the root node of the database (list)
  if (databaseRoot->contains("database")) {
    database = dynamic_cast<json::ListNode*>(&(*databaseRoot)["database"]);
  } else {
    std::cout << "DBMatDatabase: json database is ill formated (does not contain key \"database\")!"
        << std::endl;
  }
}

bool DBMatDatabase::hasDataMatrix(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  int entry_index = entryIndexByConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
        densityEstimationConfig);
  return entry_index >= 0;
}

std::string& DBMatDatabase::getDataMatrix(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  int entry_index = entryIndexByConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig);
  if (entry_index < 0) {
    // No decomposition in database
    throw sgpp::base::data_exception("Database does not contain any entry matching the "
        "decomposition");
  } else {
    json::DictNode* entry = dynamic_cast<json::DictNode*>(&((*database)[entry_index]));
    std::string& filepath = (*entry)[keyFilepath].get();
    return filepath;
  }
}

void DBMatDatabase::putDataMatrix(sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    std::string filepath, bool overwriteEntry) {
  // Check if for this setting an entry is already definied
  int entry_index = entryIndexByConfiguration(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig);
  if (entry_index < 0) {
    // Create a new list node in the json object
    json::DictNode& entry = dynamic_cast<json::DictNode&>(database->addDictValue());
    // Add a grid configuration entry
    json::DictNode& gridConfigEntry = dynamic_cast<json::DictNode&>(
        entry.addDictAttr(keyGridConfiguration));
    gridConfigEntry.addTextAttr(keyGridType, sgpp::datadriven::GeneralGridTypeParser::toString(
        gridConfig.generalType_));
    gridConfigEntry.addIDAttr(keyGridDimension, static_cast<uint64_t>(gridConfig.dim_));
    gridConfigEntry.addIDAttr(keyGridLevel, static_cast<int64_t>(gridConfig.level_));
    // Add a regularization configuration entry
    json::DictNode& regularizationConfigEntry = dynamic_cast<json::DictNode&>(entry.addDictAttr(
        keyRegularizationConfiguration));
    regularizationConfigEntry.addIDAttr(keyRegularizationStrength, regularizationConfig.lambda_);
    // Add a density estimation configuration entry
    json::DictNode& densityEstimationConfigEntry = dynamic_cast<json::DictNode&>(entry.addDictAttr(
        keyDensityEstimationConfiguration));
    densityEstimationConfigEntry.addTextAttr(keyDecompositionType,
        sgpp::datadriven::MatrixDecompositionTypeParser::toString(
            densityEstimationConfig.decomposition_));
    // Add the filepath
    entry.addTextAttr(keyFilepath, filepath);
    // Serialize the entire database
    databaseRoot->serialize(databaseFilepath);
    std::cout << "Successfully added new matrix decomposition at \"" << filepath <<
        "\" to the database" << std::endl;
  } else {
    // Update only if overwriteEntry is set to true
    if (overwriteEntry) {
      json::DictNode* entry = dynamic_cast<json::DictNode*>(&((*database)[entry_index]));
      entry->replaceTextAttr(keyFilepath, filepath);
      databaseRoot->serialize(databaseFilepath);
      std::cout << "Updated matrix decomposition to \"" << filepath << "\" in database"
          << std::endl;
    } else {
      std::cout << "Matrix decomposition already present in database. Did not update." << std::endl;
    }
  }
}

bool DBMatDatabase::gridConfigurationMatches(json::DictNode *node,
      sgpp::base::GeneralGridConfiguration& gridConfig, size_t entry_num) {
  // Check if grid general type matches
  sgpp::base::GeneralGridType gridType;
  if (node->contains(keyGridType)) {
    // Parse the grid general type
    std::string strGridType = (*node)[keyGridType].get();
    gridType = sgpp::datadriven::GeneralGridTypeParser::parse(
        strGridType);
    // Check if the general grid type matches
    if (gridConfig.generalType_ != gridType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"" << keyGridConfiguration << "\" node does not contain \"" << keyGridType <<
        "\" key and therefore is ignored!"<< std::endl;
    return false;
  }
  // Check if grid dimensionality matches
  if (node->contains(keyGridDimension)) {
    uint64_t gridDimension = (*node)[keyGridDimension].getUInt();
    if (gridConfig.dim_ != gridDimension) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"" << keyGridConfiguration << "\" node does not contain \"" << keyGridDimension <<
        "\" key and therefore is " << "ignored!" << std::endl;
    return false;
  }
  // Check if the grid level matches
  if (node->contains(keyGridLevel)) {
    // Combi grids contain a level vector of size d
    if (gridType == sgpp::base::GeneralGridType::ComponentGrid) {
      json::ListNode* entryLevelVector = dynamic_cast<json::ListNode*>(&(*node)[keyGridLevel]);
      if (entryLevelVector->size() != gridConfig.dim_) {
        std::cout << "DBMatDatabase: database entry # " << entry_num <<
            ": \"" << keyGridLevel << "\" size does not match \"" << keyGridDimension <<
            "\" key and therefore is ignored!" << std::endl;
        return false;
      }
      // Compare configuration level vector with json level vector
      sgpp::base::CombiGridConfiguration& combiGridConfig =
          dynamic_cast<sgpp::base::CombiGridConfiguration&>(gridConfig);
      if (combiGridConfig.levels.size() != gridConfig.dim_) {
        std::string what = "Invalid combi grid config: Level vector size " +
            std::to_string(combiGridConfig.levels.size()) + " does not match dimensionality of grid"
                + " configuration" + std::to_string(gridConfig.dim_);
        throw sgpp::base::data_exception(what.c_str());
      }
      for (size_t i = 0; i < gridConfig.dim_; i++) {
        json::Node* indexLevelNode = dynamic_cast<json::Node*>(&((*entryLevelVector)[i]));
        if (indexLevelNode->getInt() != combiGridConfig.levels.at(i)) {
          return false;
        }
      }
    } else {
      // Other grid types only contain a single level value
      int64_t gridLevel = (*node)[keyGridLevel].getInt();
      if (gridConfig.level_ != gridLevel) return false;
    }

  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"" << keyGridConfiguration << "\" node does not contain \"" << keyGridLevel <<
        "\" key and therefore is ignored!" << std::endl;
    return false;
  }
  // All relevant attributes of the general grid configuration match
  return true;
}

bool DBMatDatabase::regularizationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig, size_t entry_num) {
  // Check if the reguluarization strength matches
  if (node->contains(keyRegularizationStrength)) {
    double lambda = (*node)[keyRegularizationStrength].getDouble();
    if (regularizationConfig.lambda_ != lambda) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"" << keyRegularizationConfiguration << "\" node does not contain \""
        << keyRegularizationStrength << "\" key and therefore is " << "ignored!" << std::endl;
    return false;
  }
  // Regularization configuration matches
  return true;
}

bool DBMatDatabase::densityEstimationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig, size_t entry_num) {
  // Check if the decomposition type matches
  if (node->contains(keyDecompositionType)) {
    // Parse decomposition type
    std::string strDecompType = (*node)[keyDecompositionType].get();
    sgpp::datadriven::MatrixDecompositionType decompositionType =
        sgpp::datadriven::MatrixDecompositionTypeParser::parse(strDecompType);
    if (densityEstimationConfig.decomposition_ != decompositionType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num <<
        ": \"" << keyDensityEstimationConfiguration << "\" node does not contain \""
        << keyDecompositionType << "\" key and " << "therefore is ignored!" << std::endl;
    return false;
  }
  // Density estimation configuration matches
  return true;
}

int DBMatDatabase::entryIndexByConfiguration(
    sgpp::base::GeneralGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // Scan the entire database
  for (size_t i = 0; i < database->size(); i++) {
    json::DictNode* entry = dynamic_cast<json::DictNode*>(&((*database)[i]));
    // Check if the entry matches the grid configuration
    if (entry->contains(keyGridConfiguration)) {
      json::DictNode *gridConfigNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyGridConfiguration]);
      if (!gridConfigurationMatches(gridConfigNode, gridConfig, i)) continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"" << keyGridConfiguration << "\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // Check if the entry matches the regularization configuration
    if (entry->contains(keyRegularizationConfiguration)) {
      json::DictNode *regularizationConfigNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyRegularizationConfiguration]);
      if (!regularizationConfigurationMatches(regularizationConfigNode, regularizationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"" << keyRegularizationConfiguration << "\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // Check if the entry matches the density estimation configuration
    if (entry->contains(keyDensityEstimationConfiguration)) {
      json::DictNode *densityEstimationConfiNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyDensityEstimationConfiguration]);
      if (!densityEstimationConfigurationMatches(densityEstimationConfiNode,
          densityEstimationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a " <<
          "\"" << keyDensityEstimationConfiguration << "\" key and therefore is ignored!" <<
          std::endl;
      continue;
    }
    // All three configurations match -> return this entry
    if (entry->contains(keyFilepath)) {
      std::cout << "config match" << std::endl;
      return (static_cast<int>(i));
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " matches but does not contain a " <<
          "\"" << keyFilepath << "\" key and therefore is ignored!" << std::endl;
    }
  }
  return -1;
}
} /* namespace datadriven */
} /* namespace sgpp */




