// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * dataMatrixDatabase.cpp
 *
 *  Created on: Apr 22, 2018
 *      Author: dominik
 */


#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>

#include <string>

using sgpp::datadriven::DBMatDatabase;

/**
 * This example shows how to initialize a datamatrix of offline decompositions (needed for
 * online objects) which enhances the performance since the decomposition usually takes some time.
 * A database is always initialized upon a json file which contains paths to different matrix
 * decompositions identified by the configuration of the grid, adaptivity and the density
 * estimation itself.
 */


int main() {
  /**
   * Initialize the database object from a json file
   */
  std::string databasePath = "dataMatrixDatabase.json";
  DBMatDatabase database(databasePath);

  /**
   * Configure the grid
   */
  sgpp::base::CombiGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.levels.push_back(3);
  gridConfig.levels.push_back(3);

  /**
   * Configure adaptivity
   */
  sgpp::base::AdpativityConfiguration adaptivityConfig;

  /**
   * Configure regularization
   */
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 10e-7;

  /**
   * Configure density estimation
   */
  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  /**
   *  Create matrix
   */
  // std::string dbmatfilepath = "/media/d/uni/bachelor_thesis/dbmattest";
  // std::cout << "Creating dbmat" << std::endl;
  // sgpp::datadriven::DBMatOffline *db = sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
  //     gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);
  // db->buildMatrix();
  // db->decomposeMatrix();
  // db->store(dbmatfilepath);
  // std::cout << "Created dbmat" << std::endl;

  /**
   * Store the decomposition in the database. The json file is automatically updated
   */
  // database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
  //    densityEstimationConfig, dbmatfilepath, true);

  // Retrieve the file path from the database. This path can be passed to the DBMatOffline
  // factory that will instanciate a DBMatOffline object
  std::string path = database.getDataMatrix(gridConfig,
      adaptivityConfig, regularizationConfig, densityEstimationConfig);

  std::cout << "Success: " << path << std::endl;

  return 0;
}


