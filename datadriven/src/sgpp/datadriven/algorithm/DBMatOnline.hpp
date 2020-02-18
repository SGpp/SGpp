// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <list>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Class for objects that can be used in the online step of the classification
 * (The classification is divided into two parts: the offline step that does not
 * depend on the actual data and the online step that depends on the data)
 */
class DBMatOnline {
 public:
  /**
   * Constructor
   *
   * @param o a offline object
   */
  explicit DBMatOnline(DBMatOffline& o);

  DBMatOnline(const DBMatOnline& rhs) = delete;
  DBMatOnline(DBMatOnline&& rhs) = default;

  DBMatOnline& operator=(const DBMatOnline& rhs) = delete;
  DBMatOnline& operator=(DBMatOnline&& rhs) = delete;
  /**
   * Destructor
   */
  virtual ~DBMatOnline() = default;

  //  /**
  //   * Changes the weighting factor for the regularization term,
  //   * if possible (might depend on the kind of decomposition for classification)
  //   */
  void setLambda(double lambda);

  /**
   * Returns a reference to the offline object
   * @return reference to the stored offline object
   */
  DBMatOffline& getOfflineObject();

  const DBMatOffline& getOfflineObject() const;

  /**
   * Update the system matrix decomposition after the grid has been modified.
   * @param densityEstimationConfig configuration of the density estimation
   * @param grid the underlying grid
   * @param numAddedGridPoints Number of grid points inserted at the end of the grid storage
   * @param deletedGridPointIndices Indices of grid points that were deleted
   * @param lambda The last best lambda value
   * @return list of grid points, that cannot be coarsened
   */
  virtual std::vector<size_t> updateSystemMatrixDecomposition(
      DensityEstimationConfiguration& densityEstimationConfig, Grid& grid,
      size_t numAddedGridPoints, std::vector<size_t>& deletedGridPointIndices, double lambda);

 protected:
  DBMatOffline& offlineObject;
};

}  // namespace datadriven
}  // namespace sgpp
