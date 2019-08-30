/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * OperationMultipleEvalModLinearDistributed.hpp
 *
 * Created on: Mar 23, 2019
 *     Author: Jan Schopohl
 */

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalDistributed.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class is a distributed version of sgpp::base::OperationMultipleEvalModLinear, it uses
 * ScaLAPACK.
 *
 * It implements distributed OperationB for a grid with modified linear ansatzfunctions
 */
class OperationMultipleEvalModLinearDistributed : public OperationMultipleEvalDistributed {
 public:
  /**
   * Constructor
   * @param grid
   * @param dataset
   */
  OperationMultipleEvalModLinearDistributed(sgpp::base::Grid& grid, sgpp::base::DataMatrix& dataset)
      : OperationMultipleEvalDistributed(grid, dataset), storage(grid.getStorage()) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalModLinearDistributed() override {}

  /**
   * Distributed version of mult.
   * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
   *
   * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  void multDistributed(sgpp::base::DataVector& alpha, DataVectorDistributed& result) override;

  /**
   * Distributed version of multTransposed.
   * Multiplication of @f$B@f$ with vector @f$\alpha@f$
   *
   * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
   * @param result the result vector of the matrix vector multiplication
   */
  void multTransposeDistributed(
      sgpp::base::DataVector& source, DataVectorDistributed& result) override;

  double getDuration() override;

 protected:
  // reference to the GridStorage object of the grid
  sgpp::base::GridStorage& storage;
};
}  // namespace datadriven
}  // namespace sgpp
