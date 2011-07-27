/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONQUADRATUREMC_HPP
#define OPERATIONQUADRATUREMC_HPP

#include "operation/common/OperationQuadrature.hpp"
#include "grid/Grid.hpp"
#include <cstdlib>
#include <ctime>

namespace sg
{
namespace base
{

/**
 * Quadrature on any sparse grid (that has OperationMultipleEval implemented)
 * using Monte Carlo.
 * 
 */
class OperationQuadratureMC : public OperationQuadrature
{
public:
  /**
   * Constructor of OperationQuadratureMC, specifying a grid
   * object and the number of samples to use.
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param mcPaths Number of Monte Carlo samples
   *
   * @todo (pflueged) extend to error computation / arbitrary functions
   */
  OperationQuadratureMC(Grid &grid, int mcPaths);

  virtual ~OperationQuadratureMC() {}

  /**
   * Quadrature using simple MC.
   *
   * @param alpha Coefficient vector for current grid
   */
  virtual double doQuadrature(DataVector& alpha);

protected:
  // Pointer to the grid object
  Grid* grid;
  // Number of MC paths
  int mcPaths;
};

}
}

#endif /* OPERATIONQUADRATUREMC_HPP */
