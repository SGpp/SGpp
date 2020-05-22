// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

#include <iostream>
#include <string>

using sgpp::datadriven::DBMatDatabase;

/**
 * This example shows how to initialize a data matrix of offline decompositions (needed for
 * online objects) which enhances the performance since the decomposition usually takes some time.
 * A database is always initialized upon a json file which contains paths to different matrix
 * decompositions identified by the configuration of the grid, adaptivity and the density
 * estimation itself.
 *
 * NOTE (Sebastian Kreisel):
 * Use this example with a config that includes a database attribute in the
 * fitter config as well as an explicit dim in the gridConfig.
 * You may for example use config_databaseExample.json.
 *
 * From datadriven/examplesPipeline:
 * $ ./ExampleDatabase config_databaseExample.json out.txt
 *
 * writes the database to out.txt.
 */

int main(int argc, char** argv) {
  /**
   * Get the path for the configuration file and the path for the output of the decomposed matrix.
   */
  std::string configPath = "";
  std::string outputPath = "";
  if (argc != 3) {
    std::cout << "No or bad paths given, aborting" << std::endl;
    std::cout << "Usage: " << argv[0] << " pathToConfigFile pathForOutputFile" << std::endl;
    exit(1);
  } else {
    configPath = std::string{argv[1]};
    outputPath = std::string{argv[2]};
  }

  /**
   * To retrieve an offline decomposition one must define the setup the offline matrix must match.
   * This is done by defining how the sparse grid is structured, how the sparse grid performs
   * refinement and coarsening (adaptivity). Also regularization and density estimation parameters
   * must be defined to identify a decomposition. This is done by parsing a configuration file which
   * should include all these settings.
   * (gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)
   */
  sgpp::datadriven::DataMiningConfigParser parser(configPath);

  sgpp::datadriven::FitterConfigurationDensityEstimation config;
  config.readParams(parser);

  auto& gridConfig = config.getGridConfig();
  auto& adaptivityConfig = config.getRefinementConfig();
  auto& regularizationConfig = config.getRegularizationConfig();
  auto& densityEstimationConfig = config.getDensityEstimationConfig();
  auto& databaseConfig = config.getDatabaseConfig();
  /**
   * The database has to be initialized. This is done by passing the file to the json
   * database file to the constructor of the DBMatDatabase class.
   */
  DBMatDatabase database(databaseConfig.filePath_);

  /**
   * Before the matrix can be initialized the underlying grid needs to be created
   */
  std::unique_ptr<sgpp::base::Grid> grid;
  if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
    grid =
        std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createModLinearGrid(gridConfig.dim_)};
  } else if (gridConfig.type_ == sgpp::base::GridType::Linear) {
    grid = std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(gridConfig.dim_)};
  } else {
    throw sgpp::base::algorithm_exception(
        "LearnerBase::InitializeGrid: An unsupported grid type "
        "was chosen!");
  }

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.regular(gridConfig.level_);

  /**
   *  This section shows how to store a decomposition in the database. First however the
   *  matrix has to be created and decomposed. This is done using the configuration structures
   *  that the database needs to identify a decomposition.
   */
  std::cout << "Creating dbmat" << std::endl;
  sgpp::datadriven::DBMatOffline* db;
  try {
    db = sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
        gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);
  } catch (sgpp::base::factory_exception& exc)  {
    if (std::string(exc.what()).find("built without GSL") != std::string::npos) {
      std::cout << "Exception: " << exc.what() << std::endl;
      std::cout << "Skipping example..." << std::endl;
      return 0;
    } else {
      throw;
    }
  }
  db->buildMatrix(grid.get(), regularizationConfig);
  db->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  db->store(outputPath);
  std::cout << "Created dbmat" << std::endl;

  /**
   * After decomposing the data matrix it can be stored in the database. By passing the
   * configuration structures used to create the matrix the database can identify the
   * decomposition. The last parameter of the putDataMatrix method specified whether an entry with
   * the same configuration will be replaced with a new file path. Note that the database only
   * works on file paths, i.e. strings.
   */
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
                         densityEstimationConfig, outputPath, true);

  return 0;
}
