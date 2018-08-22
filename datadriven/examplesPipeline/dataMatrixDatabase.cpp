// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <string>

using sgpp::datadriven::DBMatDatabase;

/**
 * This example shows how to initialize a data matrix of offline decompositions (needed for
 * online objects) which enhances the performance since the decomposition usually takes some time.
 * A database is always initialized upon a json file which contains paths to different matrix
 * decompositions identified by the configuration of the grid, adaptivity and the density
 * estimation itself.
 */


int main() {
  /**
   * First the database has to be initialized. This is done by passing the file to the json
   * database file to the constructor of the DBMatDatabase class.
   */
  std::string databasePath = "dataMatrixDatabase.json";
  DBMatDatabase database(databasePath);

  /**
   * To retrieve an offline decomposition one must define the setup the offline matrix must match.
   * This is done by defining how the sparse grid is structured, how the sparse grid performs
   * refinement and coarsening (adaptivity). Also regularization and density estimation parameters
   * must be defined to identify a decomposition. This is done by initializing respective
   * structures for all of those settings.
   */

  /**
   * First the sparse grid itself will be specified using a structure that inherits from the
   * sgpp::base::GeneralGridConfiguration supertype. Each structure defines an unique
   * type of a sparse grid, such as a regular sparse grid, geometry aware sparse grids.
   */
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 2;
  gridConfig.level_ = 3;

  /**
   * Next the refinement and coarsening behaviour of the sparse grid is defined. In this example
   * the default values are used.
   */
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  /**
   * Also the regularization must be specified, namely the type of the regularization operator
   * and also the regularization strength, lambda. In this case the identity regularization operator
   * is chosen.
   */
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 10e-7;

  /**
   * Lastly the density estimation it self is configured. Since online / offline learning
   * will be performed the density estimation type will be
   * sgpp::datadriven::DensityEstimationType::Decomposition. Also the method of decomposition has
   * to be specified, in this case Cholesky decomposition (which supports adaptivity of the grid)
   * is chosen.
   */
  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

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
    throw sgpp::base::algorithm_exception("LearnerBase::InitializeGrid: An unsupported grid type "
        "was chosen!");
  }

  /**
   *  This section shows how to store a decomposition in the database. First however the
   *  matrix has to be created and decomposed. This is done using the configuration structures
   *  that the database needs to identify a decomposition.
   */
  std::string dbmatfilepath = "my/path";
  std::cout << "Creating dbmat" << std::endl;
  sgpp::datadriven::DBMatOffline *db = sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
      gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig);
  db->buildMatrix(grid.get(), regularizationConfig);
  db->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  db->store(dbmatfilepath);
  std::cout << "Created dbmat" << std::endl;

  /**
   * Afte decomposing the data matrix it can be stored in the database. By passing the
   * configuration structures used to create the matrix the database can identify the
   * decomposition. The last parameter of the putDataMatrix method specified whether an entry with
   * the same configuration will be replaced with a new file path. Note that the database only
   * works on file paths, i.e. strings.
   */
  database.putDataMatrix(gridConfig, adaptivityConfig, regularizationConfig,
      densityEstimationConfig, dbmatfilepath, true);

  /**
   * Lastly it is shown how to retrieve a file path from the database. To identify a decomposition
   * the entire configuration which would be used to initialize the data matrix must be passed
   * to the database. If a matching entry is found the getDataMatrix method will return the
   * file path associated with this configuration.
   */
  std::string path = database.getDataMatrix(gridConfig,
      adaptivityConfig, regularizationConfig, densityEstimationConfig);

  std::cout << "Success: " << path << std::endl;

  return 0;
}


