// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

#include <list>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::Grid;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

/**
 * Class that is used to decompose and store the left-hand-side
 * matrix for the density based classification approach
 * (The classification is divided into two parts: the offline step that does not
 * depend on the actual data and the online step that depends on the data).
 * Uses Gnu Scientific Library (GSL).
 */

class DBMatOffline {
 public:
  /**
   * Constructor
   * Build DBMatOffline Object from configuration
   *
   * @param gridConfig The configuration of the grid
   * @param adaptivityConfig The configuration of the grid adaptivity
   * @param regularizationConfig The configuration of the grid regularization
   * @param densityEstimationConfig The configuration of the matrix decomposition
   */
  explicit DBMatOffline(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::base::AdpativityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

  /**
   * Constructor
   * Create offline object from serialized offline object
   *
   * @param fileName path to the file that stores serialized offline object
   */
  explicit DBMatOffline(const std::string& fileName);

  /**
   * Copy Constructor
   *
   * @param rhs Object to copy
   */
  DBMatOffline(const DBMatOffline& rhs);

  /**
   * Default move constructor
   */
  DBMatOffline(DBMatOffline&& rhs) = default;

  /**
   * Default virtual destructor
   */
  virtual ~DBMatOffline() = default;

  /**
   * Default copy assign operator
   */
  DBMatOffline& operator=(const DBMatOffline& rhs);

  /**
   * Default move assign operator
   */
  DBMatOffline& operator=(DBMatOffline&& rhs) = default;

  /**
   * Interface for the clone idiom
   * @return a copy of this very object as a pointer to a new DBMatOffline object which is owned by
   * the caller.
   */
  virtual DBMatOffline* clone() = 0;

  /**
   * Only Offline objects based on Cholesky decomposition, or orthogonal adaptivity can be refined
   * @return true if object can be refined, else false;
   */
  virtual bool isRefineable() = 0;

  /**
   * Get a reference to the grid configuration object
   * @return Grid configuration object
   */
  sgpp::base::GeneralGridConfiguration& getGridConfig();

  /**
   * Get a reference to the grid adaptivity configuration object
   * @return Grid adaptivity configuration object
   */
  sgpp::base::AdpativityConfiguration& getAdaptivityConfig();

  /**
   * Get a reference to the grid regularization configuration object
   * @return Grid regularization configuration object
   */
  sgpp::datadriven::RegularizationConfiguration& getRegularizationConfig();

  /**
   * Get a reference to the matrix decomposition configuration object
   * @return Matrix decomposition configuration object
   */
  sgpp::datadriven::DensityEstimationConfiguration& getDensityEstimationConfig();

  /**
   * Get a reference to the decomposed matrix. Throws if matrix has not yet been decomposed.
   *
   * @return decomposed matrix
   */
  DataMatrix& getDecomposedMatrix();

  /**
   * Allows access to lhs matrix, which is meant ONLY FOR TESTING
   */
  DataMatrix& getLhsMatrix_ONLY_FOR_TESTING();

  /**
   * Returns a reference to the sparse grid
   * @return grid
   */
  Grid& getGrid();

  /**
   * Builds the right hand side matrix with or without the regularization term depending on the type
   * of decomposition
   */
  virtual void buildMatrix();

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  virtual void decomposeMatrix() = 0;

  /**
   * Prints the matrix onto standard output
   */
  void printMatrix();

  /**
   * Serialize the DBMatOffline Object
   * @param fileName path where to store the file.
   */
  virtual void store(const std::string& fileName);

  /**
   * Sets interaction Term
   * @param interactions Interaction terms used for geometrically aware grids
   */
  void setInter(std::vector<std::vector <size_t>> interactions);


 protected:
  DBMatOffline();

  // configuration of the grid
  sgpp::base::GeneralGridConfiguration gridConfig;

  // config of the grid adaptivity
  sgpp::base::AdpativityConfiguration adaptivityConfig;

  // config of the grid regularization
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;

  // config of the matrix decomposition
  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  DataMatrix lhsMatrix;              // stores the (decomposed) matrix
  bool isConstructed;                // If the matrix was built
  bool isDecomposed;                 // If the matrix was decomposed

  /**
   * An offline object works on a hierarchical basis grid.
   */
  std::unique_ptr<Grid> grid;


 public:
  // vector of interactions (if size() == 0: a regular SG is created)
  std::vector<std::vector <size_t>> interactions;

 protected:
  /**
   * Build the initial sparse grid
   */
  void InitializeGrid();

  /**
   * Read the configuration from a serialized DBMatOffline object.
   * 
   * @param fileName path of the serialized DBMatOffline object
   * @param gridConfig The configuration of the grid of the file to populate
   * @param adaptivityConfig The configuration of the grid adaptivity of the file to populate
   * @param regularizationConfig The configuration of the grid regularization of the file to populate
   * @param densityEstimationConfig The configuration of the matrix decomposition of the file to populate

   */
  void parseConfig(const std::string& fileName,
                   sgpp::base::GeneralGridConfiguration& gridConfig,
                   sgpp::base::AdpativityConfiguration& adaptivityConfig,
                   sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                   sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig) const;


  /**
   * Read the Interactionsterms from a serialized DBMatOfflibe object.
   * @param fileName path of the serialized DBMatOffline object
   * @param interactions the interactions to populate
   */
  void parseInter(const std::string& fileName,
    std::vector<std::vector<size_t>>& interactions) const;
};

}  // namespace datadriven
}  // namespace sgpp
