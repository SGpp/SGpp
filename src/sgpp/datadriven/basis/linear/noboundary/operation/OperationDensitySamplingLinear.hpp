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
      void sampling_on_all_dims(base::Grid* grid, base::DataVector* alpha, unsigned int dim_start, base::DataVector* &sampleVec);
      void sampling_on_lower_dims(base::Grid* g_in, base::DataVector* al_in, unsigned int dim_x, base::DataVector* &sampleVec, unsigned int &count);
      void sampling_on_higher_dims(base::Grid* g_in, base::DataVector* al_in, unsigned int dim_x, base::DataVector* &sampleVec, unsigned int &count);
  };

}
}
#endif /* OPERATIONDENSITYSAMPLINGLINEAR_HPP */






