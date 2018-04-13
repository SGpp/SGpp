// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYMARGINALIZEKDE_HPP
#define OPERATIONDENSITYMARGINALIZEKDE_HPP

#include <sgpp/datadriven/application/GaussianKDE.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Marginalize Probability Density Function
 */

class OperationDensityMarginalizeKDE {
 public:
  explicit OperationDensityMarginalizeKDE(datadriven::GaussianKDE& kde);
  virtual ~OperationDensityMarginalizeKDE();

  /**
   * Marginalizes (Density) functions in dimension mdim
   *
   * @param mdim marginalize in dimension mdim
   * @param marginalizedKDE marginalized kernel density
   */
  void doMarginalize(size_t mdim, datadriven::GaussianKDE& marginalizedKDE);

  /**
   * Marginalizes (Density) functions in all dimensions mdims
   * @param mdims marginalized in all dimensions in mdims
   * @param marginalizedKDE marginalized kernel density
   */
  void doMarginalize(std::vector<size_t>& mdims, datadriven::GaussianKDE& marginalizedKDE);

  /**
   * Keep applying marginalizes to (Density) Functions, until it's reduced to 1 dimension (dim_x)
   *
   * @param mdim Target dimension, all other dimensions will be marginalized
   * @param marginalizedKDE result of marginalization
   */
  void margToDimX(size_t mdim, datadriven::GaussianKDE& marginalizedKDE);

  /**
   * Keep applying marginalizes to (Density) Functions, until it's reduced to
   * the dimensions in mdims
   *
   * @param mdims Target dimensions, all other dimensions will be marginalized
   * @param marginalizedKDE result of marginalization
   */

  void margToDimXs(std::vector<size_t>& mdims, datadriven::GaussianKDE& marginalizedKDE);

 private:
  std::shared_ptr<datadriven::GaussianKDE> kde;
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONDENSITYMARGINALIZEKDE_HPP */
