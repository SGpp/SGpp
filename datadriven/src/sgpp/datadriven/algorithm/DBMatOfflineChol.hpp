/*
 * DBMatOfflineChol.hpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOfflineChol : public DBMatOfflineGE {
 public:
  DBMatOfflineChol(DBMatDensityConfiguration& oc);

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  virtual void decomposeMatrix() override;

  /**
   * Updates offline cholesky factorization based on coarsed (deletedPoints)
   * and refined (newPoints) gridPoints
   *
   * @param deletedPoints list of indices of last coarsed points
   * @param newPoints amount of refined points
   */
  void choleskyModification(size_t newPoints, std::list<size_t> deletedPoints, double lambda);

  /**
   * Updates the cholesky factor when a new grid point is added (e.g. refine)
   *
   * @param newCol DataVector with column to add to the system matrix
   * @param size columns/rows of current Cholesky factor, necessary since the
            allocated memory is increased before the Cholesky factor is modified
   */
  void choleskyAddPoint(DataVector* newCol, size_t size);

  /**
   * Permutes the rows of the cholesky factor based on permutations
   * of the system matrix (e.g. coarsening)
   *
   * @param k "left" column to permutate
   * @param l "right" column to permutate
   * @param job = 2        => left circular shift
   *	  1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ...,
   *l-1,l,k,l+1,..., size
   * 	  job = 1       => right circular shift
   * 	  1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ...,
   *l-1,l+1,...size
   */
  void choleskyPermutation(size_t k, size_t l, size_t job);
};

} /* namespace datadriven */
} /* namespace sgpp */
