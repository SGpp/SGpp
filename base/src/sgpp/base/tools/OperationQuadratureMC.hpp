// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMC_HPP
#define OPERATIONQUADRATUREMC_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Typedef for general functions that can be passed to integration methods. Requires three parameters. First, the dimensionality, then dim-many coordinates, and then further client data for the function at hand.
 */
typedef double (*FUNC)(int, double*, void*);


/**
 * Quadrature on any sparse grid (that has OperationMultipleEval implemented)
 * using Monte Carlo.
 *
 */
class OperationQuadratureMC : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureMC, specifying a grid
   * object and the number of samples to use.
   *
   * @param grid Reference to the grid object
   * @param mcPaths Number of Monte Carlo samples
   *
   */
  OperationQuadratureMC(sgpp::base::Grid& grid, int mcPaths);

  ~OperationQuadratureMC() override {}

  /**
   * Quadrature using simple MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(sgpp::base::DataVector& alpha) override;

  /**
   * Quadrature of an arbitrary function using
   * simple MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param func The function to integrate
   * @param clientdata Optional data to pass to FUNC
   */
  double doQuadratureFunc(FUNC func, void* clientdata);

  /**
   * Quadrature of the @f$L^2@f$-norm of the error,
   * @f$ ||f(x)-u(x)||_{L^2} @f$, between a given function and the
   * current sparse grid function using
   * simple MC in @f$\Omega=[0,1]^d@f$.
   *
   * @param func The function @f$f(x)@f$
   * @param clientdata Optional data to pass to FUNC
   * @param alpha Coefficient vector for current grid
   */
  double doQuadratureL2Error(FUNC func, void* clientdata,
                              sgpp::base::DataVector& alpha);

 protected:
  // Pointer to the grid object
  sgpp::base::Grid* grid;
  // Number of MC paths
  size_t mcPaths;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATUREMC_HPP */
