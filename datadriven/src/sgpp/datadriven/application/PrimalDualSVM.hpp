// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Implementation of a support vector machine in primal
 * formulation which additionally stores support vectors.
 * For non-linear classification, sparse grid kernels are applied.
 */

class PrimalDualSVM {
 public:
  /**
   * Constructor.
   *
   * @param dim The dimension of the feature space (i.e. grid size)
   * @param inputDim The dimension of the data
   * @param budget The max number of support vectors
   * @param useBias Indicates whether bias should be used
   */
  PrimalDualSVM(size_t dim, size_t inputDim, size_t budget, bool useBias);

  /**
   * Destructor.
   */
  ~PrimalDualSVM();

  /**
   * Raw prediction for a given data point and grid.
   *
   * @param grid The sparse grid which defines the transformation
   * @param x The data point
   * @param dataDim Dimension of x
   * @param trans Indicates whether x is already transformed into feature space
   * @return The raw prediction value
   */
  double predictRaw(sgpp::base::Grid& grid, sgpp::base::DataVector& x,
                    size_t dataDim, bool trans = false);

  /**
   * Class prediction for a given data point and grid.
   *
   * @param grid The sparse grid which defines the transformation
   * @param x The data point
   * @param dataDim Dimension of x
   * @return The predicted class label (-1 or 1)
   */
  int predict(sgpp::base::Grid& grid, sgpp::base::DataVector& x,
              size_t dataDim);

  /**
   * Adds a data point to the set of support vectors.
   *
   * @param grid The sparse grid which defines the transformation
   * @param x The data point
   * @param alpha The corresponding weight
   * @param dataDim Dimension of x
   */
  void add(sgpp::base::Grid& grid, sgpp::base::DataVector& x, double alpha,
           size_t dataDim);

  /**
   * Scales the normal vector w.
   *
   * @param scalar The scaling factor
   */
  void multiply(double scalar);

  // void remove(size_t idx);

  // the set of support vectors
  base::DataMatrix svs;
  // alphas corresponding to support vectors (sv-weights, signed)
  base::DataVector alphas;
  // stores the norm of each support vector in the feature space
  base::DataVector norms;
  // the normal vector (defines decision hyperplane)
  base::DataVector w;
  // normal vector computed with absolute alpha values
  base::DataVector w2;

 protected:
  // number of max support vectors to be stored
  size_t budget;
  // specifies whether bias should be applied or not
  bool useBias;
  // parameter to change position of decision hyperplane
  double bias;
};

}  // namespace datadriven
}  // namespace sgpp
