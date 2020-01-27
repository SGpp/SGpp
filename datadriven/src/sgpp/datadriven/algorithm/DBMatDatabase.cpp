// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>
#include <sgpp/datadriven/datamining/configuration/GeneralGridTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GeometryConfigurationParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GridTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

// Globally define all keys in a database entry
const std::string keyGridConfiguration = "gridConfiguration";
const std::string keyGeometryConfiguration = "geometryConfiguration";
const std::string keyRegularizationConfiguration = "regularizationConfiguration";
const std::string keyDensityEstimationConfiguration = "densityEstimationConfiguration";
const std::string keyGridType = "gridType";
const std::string keyBasisType = "basisType";
const std::string keyGridDimension = "dimension";
const std::string keyGridLevel = "level";
const std::string keyGeoDim = "dim";
const std::string keyGeoStencils = "stencils";
const std::string keyGeoStencil = "stencil";
const std::string keyGeoApply = "applyOnLayers";
const std::string keyGeoColor = "colorIndex";
const std::string keyGeoBlockLength = "blockLength";
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
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  int entry_index = entryIndexByConfiguration(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                              densityEstimationConfig);
  return entry_index >= 0;
}

bool DBMatDatabase::hasBaseDataMatrix(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  // Call with base matrix flag
  return entryIndexByConfiguration(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                   densityEstimationConfig, true) >= 0;
}

std::string& DBMatDatabase::getDataMatrix(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) {
  int entry_index = entryIndexByConfiguration(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                              densityEstimationConfig);
  if (entry_index < 0) {
    // No decomposition in database
    throw sgpp::base::data_exception(
        "Database does not contain any entry matching the "
        "decomposition");
  } else {
    json::DictNode* entry = dynamic_cast<json::DictNode*>(&((*database)[entry_index]));
    std::string& filepath = (*entry)[keyFilepath].get();
    return filepath;
  }
}

std::string& DBMatDatabase::getBaseDataMatrix(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    sgpp::base::GeneralGridConfiguration& baseGridConfig) {
  // Get index of suitable base object if existing
  int entry_index = entryIndexByConfiguration(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                              densityEstimationConfig, true);
  // If none is found throw exception
  if (entry_index < 0) {
    // No decomposition in database
    throw sgpp::base::data_exception(
        "Database does not contain any entry that can be permutated to match the "
        "decomposition");
  } else {
    // Get grid config of database entry
    json::DictNode& entry = dynamic_cast<json::DictNode&>((*database)[entry_index]);
    json::DictNode& gridConfigNode = dynamic_cast<json::DictNode&>(entry[keyGridConfiguration]);
    // Get grid dimension
    baseGridConfig.dim_ = gridConfigNode[keyGridDimension].getUInt();
    json::ListNode& entryLevelVector = dynamic_cast<json::ListNode&>(gridConfigNode[keyGridLevel]);
    // Intialize level vector
    baseGridConfig.levelVector_.clear();
    // Fill out level vector
    for (size_t i = 0; i < baseGridConfig.dim_; i++) {
      json::Node& levelNode = dynamic_cast<json::Node&>(entryLevelVector[i]);
      baseGridConfig.levelVector_.push_back(levelNode.getUInt());
    }
    // Return file path of serialized offline object
    std::string& filepath = entry[keyFilepath].get();
    return filepath;
  }
}

void DBMatDatabase::putDataMatrix(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    const std::string filepath, bool overwriteEntry) {
  // Check if for this setting an entry is already definied
  int entry_index = entryIndexByConfiguration(gridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                                              densityEstimationConfig);
  if (entry_index < 0) {
    // Create a new list node in the json object
    json::DictNode& entry = dynamic_cast<json::DictNode&>(database->addDictValue());
    // Add a grid configuration entry
    json::DictNode& gridConfigEntry =
        dynamic_cast<json::DictNode&>(entry.addDictAttr(keyGridConfiguration));
    gridConfigEntry.addTextAttr(
        keyGridType, sgpp::datadriven::GeneralGridTypeParser::toString(gridConfig.generalType_));
    gridConfigEntry.addTextAttr(
        keyBasisType, sgpp::datadriven::GridTypeParser::toString(gridConfig.type_));
    gridConfigEntry.addIDAttr(keyGridDimension, static_cast<uint64_t>(gridConfig.dim_));
    if (gridConfig.generalType_ == sgpp::base::GeneralGridType::ComponentGrid) {
      json::ListNode& level =
          dynamic_cast<json::ListNode&>(gridConfigEntry.addListAttr(keyGridLevel));
      for (size_t i = 0; i < gridConfig.levelVector_.size(); i++) {
        level.addIdValue(gridConfig.levelVector_[i]);
      }
    } else {
      gridConfigEntry.addIDAttr(keyGridLevel, static_cast<int64_t>(gridConfig.level_));
    }
    // Add a geometry configuration entry
    json::DictNode& geometryConfigEntry = 
    dynamic_cast<json::DictNode&>(entry.addDictAttr(keyGeometryConfiguration));
    json::ListNode& dims = 
    dynamic_cast<json::ListNode&>(geometryConfigEntry.addListAttr(keyGeoDim));
    for (const std::vector<int64_t>& dimension : geometryConfig.dim) {
        json::ListNode& thisDims = dynamic_cast<json::ListNode&>(dims.addListValue());
        for (const int64_t& thisDim: dimension) {
            thisDims.addIdValue(thisDim);
        }
    }
    json::ListNode& stencilsEntry = 
    dynamic_cast<json::ListNode&>(geometryConfigEntry.addListAttr(keyGeoStencils));
    for (sgpp::datadriven::StencilConfiguration stencil : geometryConfig.stencils) {
        json::DictNode& stencilEntry = 
        dynamic_cast<json::DictNode&>(stencilsEntry.addDictValue());
        stencilEntry.addTextAttr(keyGeoStencil, sgpp::datadriven::GeometryConfigurationParser::toString(stencil.stencilType));
        stencilEntry.addIDAttr(keyGeoBlockLength, static_cast<int64_t>(stencil.blockLenght));
        stencilEntry.addIDAttr(keyGeoColor, static_cast<int64_t>(stencil.colorIndex));
        json::ListNode& layersEntry =
        dynamic_cast<json::ListNode&>(stencilEntry.addListAttr(keyGeoApply));
        for (size_t i : stencil.applyOnLayers) {
            layersEntry.addIdValue(i);
        }
    }
    // Add a regularization configuration entry
    json::DictNode& regularizationConfigEntry =
        dynamic_cast<json::DictNode&>(entry.addDictAttr(keyRegularizationConfiguration));
    regularizationConfigEntry.addIDAttr(keyRegularizationStrength, regularizationConfig.lambda_);
    // Add a density estimation configuration entry
    json::DictNode& densityEstimationConfigEntry =
        dynamic_cast<json::DictNode&>(entry.addDictAttr(keyDensityEstimationConfiguration));
    densityEstimationConfigEntry.addTextAttr(
        keyDecompositionType, sgpp::datadriven::MatrixDecompositionTypeParser::toString(
                                  densityEstimationConfig.decomposition_));
    // Add the filepath
    entry.addTextAttr(keyFilepath, filepath);
    // Serialize the entire database
    databaseRoot->serialize(databaseFilepath);
    std::cout << "Successfully added new matrix decomposition at \"" << filepath
              << "\" to the database" << std::endl;
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

bool DBMatDatabase::baseGridConfigurationMatches(
    json::DictNode* node, const sgpp::base::GeneralGridConfiguration& gridConfig,
    size_t entry_num) {
  uint64_t nodeGridDimension;
  // Check if grid dimensionality matches
  if (node->contains(keyGridDimension)) {
    nodeGridDimension = (*node)[keyGridDimension].getUInt();
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridConfiguration
              << "\" node does not contain \"" << keyGridDimension << "\" key and therefore is "
              << "ignored!" << std::endl;
    return false;
  }

  json::ListNode* entryLevelVector = dynamic_cast<json::ListNode*>((&(*node)[keyGridLevel]));
  if (entryLevelVector == nullptr) {
    return false;
  }
  if (entryLevelVector->size() != nodeGridDimension) {
    return false;
  }
  std::vector<size_t> levelVec;

  for (size_t i = 0; i < nodeGridDimension; i++) {
    json::Node* levelNode = dynamic_cast<json::Node*>(&((*entryLevelVector)[i]));
    levelVec.push_back(levelNode->getUInt());
  }

  return PermutationUtil::isPermutation(
      PermutationUtil::deleteOnesFromLevelVec(levelVec),
      PermutationUtil::deleteOnesFromLevelVec(gridConfig.levelVector_));
}

bool DBMatDatabase::gridConfigurationMatches(json::DictNode* node,
                                             const sgpp::base::GeneralGridConfiguration& gridConfig,
                                             size_t entry_num) {
  // Check if grid general type matches
  sgpp::base::GeneralGridType gridType;
  if (node->contains(keyGridType)) {
    // Parse the grid general type
    std::string strGridType = (*node)[keyGridType].get();
    gridType = sgpp::datadriven::GeneralGridTypeParser::parse(strGridType);
    // Check if the general grid type matches
    if (gridConfig.generalType_ != gridType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridConfiguration
              << "\" node does not contain \"" << keyGridType << "\" key and therefore is ignored!"
              << std::endl;
    return false;
  }
  if (node->contains(keyBasisType)) {
    // Parse the basis type
    std::string strBasisType = (*node)[keyBasisType].get();
    sgpp::base::GridType basisType = sgpp::datadriven::GridTypeParser::parse(strBasisType);
    // Check if the basis type matches
    if (gridConfig.type_ != basisType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridConfiguration
              << "\" node does not contain \"" << keyBasisType << "\" key and therefore is ignored!"
              << std::endl;
    return false;
  }
  // Check if grid dimensionality matches
  if (node->contains(keyGridDimension)) {
    uint64_t gridDimension = (*node)[keyGridDimension].getUInt();
    if (gridConfig.dim_ != gridDimension) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridConfiguration
              << "\" node does not contain \"" << keyGridDimension << "\" key and therefore is "
              << "ignored!" << std::endl;
    return false;
  }
  // Check if the grid level matches
  if (node->contains(keyGridLevel)) {
    // Combi grids contain a level vector of size d
    if (gridType == sgpp::base::GeneralGridType::ComponentGrid) {
      json::ListNode* entryLevelVector = dynamic_cast<json::ListNode*>(&(*node)[keyGridLevel]);
      if (entryLevelVector->size() != gridConfig.dim_) {
        std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridLevel
                  << "\" size does not match \"" << keyGridDimension
                  << "\" key and therefore is ignored!" << std::endl;
        return false;
      }
      // Compare configuration level vector with json level vector
      if (gridConfig.levelVector_.size() != gridConfig.dim_) {
        std::string what = "Invalid combi grid config: Level vector size " +
                           std::to_string(gridConfig.levelVector_.size()) +
                           " does not match dimensionality of grid" + " configuration" +
                           std::to_string(gridConfig.dim_);
        throw sgpp::base::data_exception(what.c_str());
      }
      for (size_t i = 0; i < gridConfig.dim_; i++) {
        json::Node* indexLevelNode = dynamic_cast<json::Node*>(&((*entryLevelVector)[i]));
        if (indexLevelNode->getUInt() != gridConfig.levelVector_.at(i)) {
          return false;
        }
      }
    } else {
      // Other grid types only contain a single level value
      int64_t gridLevel = (*node)[keyGridLevel].getInt();
      if (gridConfig.level_ != gridLevel) return false;
    }

  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGridConfiguration
              << "\" node does not contain \"" << keyGridLevel << "\" key and therefore is ignored!"
              << std::endl;
    return false;
  }
  // All relevant attributes of the general grid configuration match
  return true;
}

bool DBMatDatabase::geometryConfigurationMatches(
    json::DictNode* node, const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    size_t entry_num) {
    
  // Check if dims match
  if (node->contains(keyGeoDim)) {
    json::ListNode* entryDims = dynamic_cast<json::ListNode*>(&(*node)[keyGeoDim]);
    if (entryDims->size() != geometryConfig.dim.size()) {
        std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGeoDim
                  << "\" size does not match and therefore is ignored!" << std::endl;
        return false;
    }
    for (size_t i = 0; i < entryDims->size(); i++) {
        json::ListNode* entryDim = dynamic_cast<json::ListNode*>(&((*entryDims)[i]));
        if (entryDim->size() != geometryConfig.dim[i].size()) {
            return false;
        }
        for (size_t j = 0; j < entryDim->size(); j++) {
            json::Node* indexDimNode = dynamic_cast<json::Node*>(&((*entryDim)[j]));
            if (indexDimNode->getInt() != geometryConfig.dim[i].at(j)) {
                return false;
            }
        }
    }
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGeometryConfiguration 
    << "\" node does not contain \"" << keyGeoDim << "\" key and therefore is ignored!"
    << std::endl;
    return false;
  }
  // Check if stencils match
  if (node->contains(keyGeoStencils)) {
    json::ListNode* entryStencils = dynamic_cast<json::ListNode*>(&(*node)[keyGeoStencils]);
    if (entryStencils->size() != geometryConfig.stencils.size()) {
        std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGeoStencils 
        << "\" size does not match and therefore is ignored!" << std::endl;
        return false;
    }
    for (size_t i = 0; i < entryStencils->size(); i++) {
        json::DictNode* entryStencil = dynamic_cast<json::DictNode*>(&((*entryStencils)[i]));
        if (entryStencil->contains(keyGeoColor)) {
            int64_t colorIndex = (*entryStencil)[keyGeoColor].getInt();
            if (geometryConfig.stencils[i].colorIndex != colorIndex) return false;
        } else {
            return false;
        }
        if (entryStencil->contains(keyGeoBlockLength)) {
            if (geometryConfig.stencils[i].blockLenght != (*entryStencil)[keyGeoBlockLength].getUInt()) return false;
        } else {
            return false;
        }
        if (entryStencil->contains(keyGeoStencil)) {
            // Parse stencil type
            std::string strStencilType = (*entryStencil)[keyGeoStencil].get();
            sgpp::datadriven::StencilType stencilType = 
            sgpp::datadriven::GeometryConfigurationParser::parseStencil(strStencilType);
            if (geometryConfig.stencils[i].stencilType != stencilType) return false;
        } else {
            return false;
        }
        if (entryStencil->contains(keyGeoApply)) {
            json::ListNode* entryApply = dynamic_cast<json::ListNode*>(&(*entryStencil)[keyGeoApply]);
            if (entryApply->size() != geometryConfig.stencils[i].applyOnLayers.size()) {
                return false;
            }
            for (size_t j = 0; j < entryApply->size(); j++) {
                json::Node* indexApplyNode = dynamic_cast<json::Node*>(&((*entryApply)[j]));
                if (indexApplyNode->getUInt() != geometryConfig.stencils[i].applyOnLayers.at(j)) {
                    return false;
                }
            }
        } else {
            return false;
        }
    }
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \"" << keyGeometryConfiguration 
    << "\" node does not contain \"" << keyGeoStencils << "\" key and therefore is ignored!"
    << std::endl;
    return false;
  }
  // Geometry Matches!
  return true;
}


bool DBMatDatabase::regularizationConfigurationMatches(
    json::DictNode* node, const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    size_t entry_num) {
  // Check if the reguluarization strength matches
  if (node->contains(keyRegularizationStrength)) {
    double lambda = (*node)[keyRegularizationStrength].getDouble();
    if (regularizationConfig.lambda_ != lambda) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \""
              << keyRegularizationConfiguration << "\" node does not contain \""
              << keyRegularizationStrength << "\" key and therefore is "
              << "ignored!" << std::endl;
    return false;
  }
  // Regularization configuration matches
  return true;
}

bool DBMatDatabase::densityEstimationConfigurationMatches(
    json::DictNode* node,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    size_t entry_num) {
  // Check if the decomposition type matches
  if (node->contains(keyDecompositionType)) {
    // Parse decomposition type
    std::string strDecompType = (*node)[keyDecompositionType].get();
    sgpp::datadriven::MatrixDecompositionType decompositionType =
        sgpp::datadriven::MatrixDecompositionTypeParser::parse(strDecompType);
    if (densityEstimationConfig.decomposition_ != decompositionType) return false;
  } else {
    std::cout << "DBMatDatabase: database entry # " << entry_num << ": \""
              << keyDensityEstimationConfiguration << "\" node does not contain \""
              << keyDecompositionType << "\" key and "
              << "therefore is ignored!" << std::endl;
    return false;
  }
  // Density estimation configuration matches
  return true;
}

int DBMatDatabase::entryIndexByConfiguration(
    const sgpp::base::GeneralGridConfiguration& gridConfig,
    const sgpp::datadriven::GeometryConfiguration& geometryConfig,
    const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    bool findBaseConfig) {
  if (findBaseConfig && gridConfig.generalType_ != sgpp::base::GeneralGridType::ComponentGrid) {
    throw sgpp::base::algorithm_exception(
        "Base matrices can only be found for anisotrophic grids.");
  }
  // Scan the entire database
  for (size_t i = 0; i < database->size(); i++) {
    json::DictNode* entry = dynamic_cast<json::DictNode*>(&((*database)[i]));
    // Check if the entry matches the grid configuration
    if (entry->contains(keyGridConfiguration)) {
      json::DictNode* gridConfigNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyGridConfiguration]);
      if (findBaseConfig) {
        if (!baseGridConfigurationMatches(gridConfigNode, gridConfig, i)) continue;
      } else {
        if (!gridConfigurationMatches(gridConfigNode, gridConfig, i)) continue;
      }
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a "
                << "\"" << keyGridConfiguration << "\" key and therefore is ignored!" << std::endl;
      continue;
    }
    // Check if the entry matches the geometry configuration
    if (entry->contains(keyGeometryConfiguration)) {
      json::DictNode* geometryConfigNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyGeometryConfiguration]);
      if (!geometryConfigurationMatches(geometryConfigNode, geometryConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a "
                << "\"" << keyGeometryConfiguration << "\" key and therefore is ignored!"
                << std::endl;
      continue;
    }
    // Check if the entry matches the regularization configuration
    if (entry->contains(keyRegularizationConfiguration)) {
      json::DictNode* regularizationConfigNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyRegularizationConfiguration]);
      if (!regularizationConfigurationMatches(regularizationConfigNode, regularizationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a "
                << "\"" << keyRegularizationConfiguration << "\" key and therefore is ignored!"
                << std::endl;
      continue;
    }
    // Check if the entry matches the density estimation configuration
    if (entry->contains(keyDensityEstimationConfiguration)) {
      json::DictNode* densityEstimationConfiNode =
          dynamic_cast<json::DictNode*>(&(*entry)[keyDensityEstimationConfiguration]);
      if (!densityEstimationConfigurationMatches(densityEstimationConfiNode,
                                                 densityEstimationConfig, i))
        continue;
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " does not contain a "
                << "\"" << keyDensityEstimationConfiguration << "\" key and therefore is ignored!"
                << std::endl;
      continue;
    }
    // All three configurations match -> return this entry
    if (entry->contains(keyFilepath)) {
      std::cout << "config match" << std::endl;
      return (static_cast<int>(i));
    } else {
      std::cout << "DBMatDatabase: database entry # " << i << " matches but does not contain a "
                << "\"" << keyFilepath << "\" key and therefore is ignored!" << std::endl;
    }
  }
  return -1;
}
} /* namespace datadriven */
} /* namespace sgpp */
