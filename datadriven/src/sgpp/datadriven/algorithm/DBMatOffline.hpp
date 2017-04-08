// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMATOFFLINE_H_
#define DBMATOFFLINE_H_

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

#include <gsl/gsl_permutation.h>

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
   *
   * @param oc configuration for this offline object
   */
  explicit DBMatOffline(const DBMatDensityConfiguration& oc);

  /**
   * Constructor
   *
   * @param fname name of the file that stores the matrix + configuration
   */
  explicit DBMatOffline(const std::string& fname);

  /**
   * Copy Constructor
   *
   * @param rhs Object to copy
   */
  DBMatOffline(const DBMatOffline& rhs);

  DBMatOffline(DBMatOffline&& rhs) = default;

  virtual ~DBMatOffline() = default;

  DBMatOffline& operator=(const DBMatOffline& rhs);

  DBMatOffline& operator=(DBMatOffline&& rhs) = default;

  virtual DBMatOffline* clone() = 0;

  /**
   * Returns a pointer to the configuration
   */
  DBMatDensityConfiguration& getConfig();

  /**
   * Returns a pointer to the decomposed matrix
   */
  DataMatrix& getDecomposedMatrix();

  /**
   * Returns a reference to the sparse grid
   * (if this offline object uses a sparse grid, otherwise NULL)
   */
  Grid& getGrid();

  /**
   * Builds the right hand side matrix with or without the regularization term
   * depending
   * on the type of decomposition
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
   * Stores the matrix + configuration
   *
   * @param fname the file name
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
   * Method to initialize a sparse grid
   */
  void InitializeGrid();

  void parseConfig(const std::string& fileName, DBMatDensityConfiguration& config) const;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DBMATOFFLINE_H_ */

#endif /* USE_GSL */
