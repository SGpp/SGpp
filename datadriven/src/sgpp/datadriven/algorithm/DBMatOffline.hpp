// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/ParallelConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModifiedLinear.hpp>

#include <list>
#include <memory>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;

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
   * Get a reference to the decomposed matrix. Throws if matrix has not yet been decomposed.
   *
   * @return decomposed matrix
   */
  DataMatrix& getDecomposedMatrix();

  /**
   * Get a reference to the inverse matrix
   *
   * @return inverse matrix
   */
  DataMatrix& getInverseMatrix();

  /**
   * Get a reference to the distributed decomposed matrix. Throws if matrix has not yet been
   * decomposed. In order to return valid data, syncDistributedDecomposition() has to be called if
   * the decomposition was changed (after refinements).
   */
  DataMatrixDistributed& getDecomposedMatrixDistributed();

  /**
   * Get a reference to the distributed decomposed matrix, analogously to
   * getDecomposedMatrixDistributed()
   */
  DataMatrixDistributed& getDecomposedInverseDistributed();

  /**
   * Synchronizes the decomposed matrix.
   * Override if more matrices have to be synched.
   * @param processGrid process grid to distribute the matrix on
   * @param parallelConfig
   */
  virtual void syncDistributedDecomposition(std::shared_ptr<BlacsProcessGrid> processGrid,
                                            const ParallelConfiguration& parallelConfig);

  /**
   * Synchronizes the inverse matrix
   * @param processGrid process grid to distribute the matrix on
   * @param parallelConfig
   */
  virtual void syncDistributedInverse(std::shared_ptr<BlacsProcessGrid> processGrid,
                                      const ParallelConfiguration& parallelConfig);

  /**
   * Allows access to lhs matrix, which is meant ONLY FOR TESTING
   */
  DataMatrix& getLhsMatrix_ONLY_FOR_TESTING();

  /**
   * Builds the right hand side matrix with or without the regularization term depending on the type
   * of decomposition
   * @param grid The grid object the matrix is based on
   * @param regularizationConfig Configures the regularization which is incorporated into the lhs
   */
  virtual void buildMatrix(Grid* grid, RegularizationConfiguration& regularizationConfig);

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   * @param regularizationConfig the regularization configuration
   * @param densityEstimationConfig the density estimation configuration
   */
  virtual void decomposeMatrix(RegularizationConfiguration& regularizationConfig,
                               DensityEstimationConfiguration& densityEstimationConfig) = 0;

  /**
   * Prints the matrix onto standard output
   */
  void printMatrix();

  /**
   * computes vectors of L2 products of grid points for refinement
   * @param mat_refine matrix to store L2 vectors
   * @param grid underlying grid
   * @param newPoints amount of points to refine
   */
  void compute_L2_refine_vectors(DataMatrix* mat_refine, Grid* grid, size_t newPoints);

  /*
   * explicitly computes the inverse of the decomposed offline matrix
   * @param inv the matrix to store the computed inverse
   */
  virtual void compute_inverse();

  /**
   * Serialize the DBMatOffline Object
   * @param fileName path where to store the file.
   */
  virtual void store(const std::string& fileName);

  /**
   * Returns the dimensionality of the quadratic lhs matrix (i.e. the number of rows)
   * @return the grid size
   */
  virtual size_t getGridSize();

  /**
   * Returns the decomposition type of the DBMatOffline object
   * @return the type of matrix decomposition
   */
  virtual sgpp::datadriven::MatrixDecompositionType getDecompositionType() = 0;

 protected:
  DBMatOffline();
  DataMatrix lhsMatrix;   // stores the (decomposed) matrix
  bool isConstructed;     // If the matrix was built
  bool isDecomposed;      // If the matrix was decomposed
  DataMatrix lhsInverse;  // stores the explicitly computed inverse (only in SMW case)

  // distributed lhs, only initialized in ScaLAPACK version
  DataMatrixDistributed lhsDistributed;
  DataMatrixDistributed lhsDistributedInverse;

 public:
  // vector of interactions (if size() == 0: a regular SG is created)
  std::vector<std::vector<size_t>> interactions;

 protected:
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
