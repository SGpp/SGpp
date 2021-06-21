// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <iostream>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class provides the quadrature of a sparse grid function
 */
class OperationWeightedQuadrature {
 public:
  /**
   * Constructor
   */
  OperationWeightedQuadrature() {}

  /**
   * Destructor
   */
  virtual ~OperationWeightedQuadrature() {}

  /**
   *Integrate the sparse grid function w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdfs			probability density functions
   */
  virtual double doWeightedQuadrature(DataVector& alpha, sgpp::base::DistributionsVector pdfs) = 0;
};

}  // namespace base
}  // namespace sgpp
