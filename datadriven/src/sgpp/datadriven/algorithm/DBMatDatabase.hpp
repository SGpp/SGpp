// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * A database class to store and retrieve online matrix decompositions for the sparse grid
 * density estimation. The class works on a json file.
 */
class DBMatDatabase {
 public:
  /**
   * Initializes the database from a json filepath.
   * @param filepath the path to the json database
   */
  explicit DBMatDatabase(const std::string& filepath);
  virtual ~DBMatDatabase() = default;

  /**
   * Scans the entire database and checks whether any entry matches the configuration.
   * @param gridConfig the grid configuration the matrix must match
   * @param adaptivityConfig the adaptivity configuration the matrix must match
   * @param regularizationConfig the regularization configuration the matrix must match
   * @param densityEstimationConfig the density estimation configuration the matrix must match
   * @return weather the configuration is held in the database
   */
  bool hasDataMatrix(sgpp::base::GeneralGridConfiguration& gridConfig,
                     sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                     sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                     sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);
  /**
   * Scans the entire database and checks whether there exists an entry that can be permutated to
   * match the configurations. A entry matches if all config properties except the gridConfig's
   * level vector are equal and level vector elements unequal 1 are set equal.
   * @param gridConfig the grid configuration the matrix must match
   * @param adaptivityConfig the adaptivity configuration the matrix must match
   * @param regularizationConfig the regularization configuration the matrix must match
   * @param densityEstimationConfig the density estimation configuration the matrix must match
   * @return weather the configuration is held in the database
   */
  bool hasBaseDataMatrix(sgpp::base::GeneralGridConfiguration& gridConfig,
                         sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                         sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

  /**
   * Scans the entire database and finds the first entry that matches the configurations.
   * @param gridConfig the grid configuration the matrix must match
   * @param adaptivityConfig the adaptivity configuration the matrix must match
   * @param regularizationConfig the regularization configuration the matrix must match
   * @param densityEstimationConfig the density estimation configuration the matrix must match
   * @return Returns the string of the datamatrix if any match was obtained and throws an exception
   * otherwise
   */
  std::string& getDataMatrix(
      sgpp::base::GeneralGridConfiguration& gridConfig,
      sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

  /**
   * Scans the entire database and returns the first entry that can be permutated to match the
   * configurations. A entry matches if all config properties except the gridConfig's level vector
   * are equal and level vector elements unequal 1 are set equal.
   * @param gridConfig the grid configuration the matrix must match
   * @param adaptivityConfig the adaptivity configuration the matrix must match
   * @param regularizationConfig the regularization configuration the matrix must match
   * @param densityEstimationConfig the density estimation configuration the matrix must match
   * @param baseGridConfig is overriten with the entry's grid config 
   * @return Returns the string of the datamatrix if any match was obtained and throws an exception
   * otherwise
   */
  std::string& getBaseDataMatrix(
      sgpp::base::GeneralGridConfiguration& gridConfig,
      sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      sgpp::base::GeneralGridConfiguration baseGridConfig);

  /**
   * Puts a filepath for a given configuration in the database. The filepath refers to the matrix
   * file. If for this configuration a filepath is already present in the database the filepath
   * is updated if and only if the overwriteEntry parameter is set (default = false).
   * @param gridConfig the grid configuration the matrix matches
   * @param adaptivityConfig the adaptivity configuration the matrix matches
   * @param regularizationConfig the regularization configuration the matrix matches
   * @param densityEstimationConfig the density estimation configuration the matrix matches
   * @param filepath the path where the matrix decomposition is located at
   * @param overwriteEntry replaces existing entries with the same configuration if and only if
   * this parameter is set
   */
  void putDataMatrix(sgpp::base::GeneralGridConfiguration& gridConfig,
                     sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                     sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                     sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
                     std::string filepath, bool overwriteEntry = false);

 private:
  /**
   * Path to the json file containing the database
   */
  std::string databaseFilepath;

  /**
   * Root json list node containing database entries
   */
  json::ListNode* database;

  /**
   * Root node of the json file (since appearantly the json module does not support a list node
   * as root node of a json file)
   */
  std::unique_ptr<json::JSON> databaseRoot;

  /**
   * Scans the entire database and finds the first entry that matches the configurations. Returns
   * the index of the entry in the database ListNode or -1 if no entry matches.
   * @param gridConfig the grid configuration the matrix matches
   * @param adaptivityConfig the adaptivity configuration the matrix matches
   * @param regularizationConfig the regularization configuration the matrix matches
   * @param densityEstimationConfig the density estimation configuration the matrix matches
   * @return the index of the entry that matches the configuration or -1 if not entry matches
   */
  int entryIndexByConfiguration(
      sgpp::base::GeneralGridConfiguration& gridConfig,
      sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      bool findBaseConfig = false);

  /**
   * Checks weather the grid configuration of a json dict node representing a database entry root
   * matches the grid configuration passed to the database.
   * @param node the root node of the grid configuration
   * @param gridConfig the configuration to check
   * @param entry_num the index of the entry the root node belongs to
   * @return if the passed grid configuration matches the grid configuration of the node
   */
  bool gridConfigurationMatches(json::DictNode* node,
                                sgpp::base::GeneralGridConfiguration& gridConfig, size_t entry_num);

  bool baseGridConfigurationMatches(json::DictNode* node,
                                    sgpp::base::GeneralGridConfiguration& gridConfig,
                                    size_t entry_num);

  /**
   * Checks weather the regularization configuration of a json dict node representing a
   * database entry root matches the regularization configuration passed to the database.
   * @param node the root node of the regularization configuration
   * @param regularizationConfig the configuration to check
   * @param entry_num the index of the entry the root node belongs to
   * @return if the passed regularization configuration matches the regularization configuration
   * of the node
   */
  bool regularizationConfigurationMatches(
      json::DictNode* node, sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      size_t entry_num);

  /**
   * Checks weather the regularization configuration of a json dict node representing a
   * database entry root matches the regularization configuration passed to the database.
   * @param node the root node of the density estimation configuration
   * @param densityEstimationConfig the configuration to check
   * @param entry_num the index of the entry the root node belongs to
   * @return if the passed density estimation configuration matches the density estimation
   * configuration of the node
   */
  bool densityEstimationConfigurationMatches(
      json::DictNode* node,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig, size_t entry_num);
};
} /* namespace datadriven */
} /* namespace sgpp */
