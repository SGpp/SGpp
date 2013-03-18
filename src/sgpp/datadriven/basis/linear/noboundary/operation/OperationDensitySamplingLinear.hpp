/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#ifndef OPERATIONDENSITYSAMPLINGLINEAR_HPP
#define OPERATIONDENSITYSAMPLINGLINEAR_HPP

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensitySampling.hpp"

namespace sg
{
namespace datadriven
{

  /**
   * keep applying marginalize to function until it's reduced to only 1 dimension
   */

  class OperationDensitySamplingLinear : public OperationDensitySampling
  {
    public:
  	  OperationDensitySamplingLinear(base::Grid* grid) : grid(grid) {}
	  virtual ~OperationDensitySamplingLinear() {}

	  /**
	   * Keep applying marginalizes to (Density) Functions, until it's reduced to 1 dimension (dim_x)
	   *
	   * @param alpha Coefficient vector for current grid
	   * @param dim_x Target dimension number, all dimensions higher and lower than dim_x will be removed
	   * @param dim_h Highest dimension number of grid
	   */
	  void doSampling(base::DataVector* alpha, base::DataMatrix* &samples, size_t num_samples);

    protected:
      base::Grid* grid;
      void doSampling_start_dimX(base::Grid* g_in, base::DataVector* a_in, size_t dim_start, base::DataVector* &sampleVec);
      void doSampling_in_next_dim(base::Grid* g_in, base::DataVector* a_in, size_t dim_x, base::DataVector* &sampleVec, size_t &curr_dim);
  };

}
}
#endif /* OPERATIONDENSITYSAMPLINGLINEAR_HPP */






