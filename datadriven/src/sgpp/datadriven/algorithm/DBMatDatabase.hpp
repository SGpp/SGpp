// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * DBMatDatabase.hpp
 *
 *  Created on: Apr 22, 2018
 *      Author: dominik
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_DBMATDATABASE_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_DBMATDATABASE_HPP_

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class DBMatDatabase{
 public:
  explicit DBMatDatabase(const std::string& filepath);
  virtual ~DBMatDatabase() = default;

  /**
   * Scans the entire database and finds the first entry that matches the configurations. A matrix
   * is created based on the precomputed matrix decomposition the entry is referencing. If no match
   * is obtained then null is returned.
   */
  DBMatOffline* dataMatrixFromDatabase(sgpp::base::GeneralGridConfiguration& gridConfig,
      sgpp::base::AdpativityConfiguration& adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

 private:
  json::ListNode* database;
  std::unique_ptr<json::JSON> databaseRoot;

  /**
   * Scans the entire database and finds the first entry that matches the configurations. If
   * none matches null is returned. Note that not the "best" but the "first" match is returned.
   */
  std::string matrixFileByConfiguration(sgpp::base::GeneralGridConfiguration& gridConfig,
      sgpp::base::AdpativityConfiguration& adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

  /**
   * Checks weather the grid configuration of a json dict node representing a database entry root
   * matches the grid configuration passed to the database.
   */
  bool gridConfigurationMatches(json::DictNode *node,
      sgpp::base::GeneralGridConfiguration& gridConfig, size_t entry_num);

  /**
   * Checks weather the regularization configuration of a json dict node representing a
   * database entry root matches the regularization configuration passed to the database.
   */
  bool regularizationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::RegularizationConfiguration& gridConfig, size_t entry_num);

  /**
   * Checks weather the regularization configuration of a json dict node representing a
   * database entry root matches the regularization configuration passed to the database.
   */
  bool densityEstimationConfigurationMatches(json::DictNode *node,
      sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      size_t entry_num);
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_ALGORITHM_DBMATDATABASE_HPP_ */
