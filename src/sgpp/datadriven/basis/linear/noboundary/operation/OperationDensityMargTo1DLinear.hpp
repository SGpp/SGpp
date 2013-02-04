/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#ifndef OPERATIONDENSITYMARGTO1DLINEAR_HPP
#define OPERATIONDENSITYMARGTO1DLINEAR_HPP

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensityMargTo1D.hpp"
#include <cstring>

namespace sg
{
namespace datadriven
{

  /**
   * keep applying marginalize to function until it's reduced to only 1 dimension
   */

  class OperationDensityMargTo1DLinear : public OperationDensityMargTo1D
  {
    public:
  	  OperationDensityMargTo1DLinear(base::Grid* grid) : grid(grid) {}
	  virtual ~OperationDensityMargTo1DLinear() {}

	  /**
	   * Keep applying marginalizes to (Density) Functions, until it's reduced to 1 dimension (dim_x)
	   *
	   * @param alpha Coefficient vector for current grid
	   * @param grid_x Referenz of grid pointer
	   * @param alpha_x Coefficient vector for new grid (grid_x). Will be resized.
	   * @param dim_x Target dimension number, all dimensions higher and lower than dim_x will be removed
	   * @param dim_h Highest dimension number of grid
	   */
	  void margToDimX(base::DataVector* alpha, base::Grid* &grid_x, base::DataVector* &alpha_x, unsigned int dim_x);

    protected:
      base::Grid* grid;
      void margRemoveLowerDims(base::Grid* g_in, base::DataVector* al_in, base::Grid* &g_out, base::DataVector* &al_out, unsigned int dim_x, unsigned int dim_h);
      void margToDim0(base::Grid* g_in, base::DataVector* al_in, base::Grid* &g_out, base::DataVector* &al_out);
  };

}
}
#endif /* OPERATIONDENSITYMARGTO1DLINEAR_HPP */





