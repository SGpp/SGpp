// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIER_HPP
#define SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>

#include <vector>

namespace SGPP {
namespace optimization {

/**
 * Abstract operation for hierarchization and dehierarchization for
 * multiple sets of function values at the grid nodes.
 */
class OperationMultipleHierarchisation {
 public:
  /**
   * Constructor.
   */
  OperationMultipleHierarchisation() {}

  /**
   * Destructor.
   */
  virtual ~OperationMultipleHierarchisation() {}

  /**
   * Virtual method for hierarchizing for one set of function values.
   *
   * @param[in,out] nodeValues before: vector of function values at
   *                           the grid points,
   *                           after: vector of hierarchical coefficients
   * @return                   whether hierarchisation was successful
   */
  virtual bool doHierarchisation(base::DataVector& nodeValues) = 0;

  /**
   * Virtual method for dehierarchizing for one set of function values.
   *
   * @param[in,out] alpha before: vector of hierarchical coefficients,
   *                      after: vector of function values at
   *                      the grid points
   */
  virtual void doDehierarchisation(base::DataVector& alpha) = 0;

  /**
   * Pure virtual method for hierarchizing for multiple sets of
   * function values.
   *
   * @param[in,out] nodeValues before: matrix of function values at
   *                           the grid points,
   *                           after: matrix of hierarchical coefficients
   * @return                   whether hierarchisation was successful
   */
  virtual bool doHierarchisation(base::DataMatrix& nodeValues) = 0;

  /**
   * Pure virtual method for dehierarchizing for multiple sets of
   * coefficients.
   *
   * @param[in,out] alpha before: matrix of hierarchical coefficients,
   *                      after: matrix of function values at
   *                      the grid points
   */
  virtual void doDehierarchisation(base::DataMatrix& alpha) = 0;
};
}  // namespace optimization
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIER_HPP */
