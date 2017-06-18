// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

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
   * @param config configuration for this offline object
   */
  explicit DBMatOffline(const DBMatDensityConfiguration& config);

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
   * Only Offline objects based on Cholesky decomposition can be refined
   * @return true if object can be refined, else false;
   */
  virtual bool isRefineable() = 0;

  /**
   * Get a reference to the configuration object
   * @return Configuration object
   */
  DBMatDensityConfiguration& getConfig();

  /**
   * Get a reference to the decomposed matrix. Throws if matrix has not yet been decomposed.
   * @return decomposed matrix
   */
  DataMatrix& getDecomposedMatrix();

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

 protected:
  DBMatOffline();

  DBMatDensityConfiguration config;  // configuration for this offline object
  DataMatrix lhsMatrix;              // stores the (decomposed) matrix
  bool isConstructed;                // If the matrix was built
  bool isDecomposed;                 // If the matrix was decomposed

  /**
   * An offline object works on a hierarchical basis grid.
   */
  std::unique_ptr<Grid> grid;

  /**
   * Build the initial sparse grid
   */
  void InitializeGrid();

  /**
   * Read the DBMatDensityConfiguration object from a serialized DBMatOfflibe object.
   * @param fileName path of the serialized DBMatOffline object
   * @param config the configuration file to populate
   */
  void parseConfig(const std::string& fileName, DBMatDensityConfiguration& config) const;
};

}  // namespace datadriven
}  // namespace sgpp
