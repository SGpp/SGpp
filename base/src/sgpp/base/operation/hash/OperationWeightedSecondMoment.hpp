// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/Distribution.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class provides the second moment of a sparse grid function f w.r.t. a probability density
 * function \rho, i.e.
 * \int f(x)^2 \rho (x) dx
 */
class OperationWeightedSecondMoment {
 public:
  /**
   * Constructor
   */
  OperationWeightedSecondMoment() {}

  /**
   * Destructor
   */
  virtual ~OperationWeightedSecondMoment() {}

  /**
   * Integrate the sparse grid function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdf			probability density function
   * @parm quadOrder	order for the gauss Legendre quadrature
   */
  virtual double doWeightedQuadrature(DataVector& alpha,
                                      std::shared_ptr<sgpp::base::Distribution> pdf,
                                      size_t quadOrder) = 0;
};

}  // namespace base
}  // namespace sgpp
