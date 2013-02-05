/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#ifndef OPERATIONDENSITYSAMPLING_HPP_
#define OPERATIONDENSITYSAMPLING_HPP_

#include "base/grid/Grid.hpp"

namespace sg
{
namespace datadriven
{

  /**
   * Sampling on all dimensions
   */

  class OperationDensitySampling
  {
    public:
      OperationDensitySampling() {}
      virtual ~OperationDensitySampling() {}

	  /**
	   * Draw samples on all dimensions.
	   * Each dimension x as a base, a number of samples would be draw based on dim_x
	   *
	   * @param alpha Coefficient vector for current grid
	   * @param samples output matrix
	   * @param num_sampels number of samples to draw
	   */
      virtual void doSampling(base::DataVector* alpha, base::DataMatrix* &samples, size_t num_samples) = 0;
  };

}
}
#endif /* OPERATIONDENSITYSAMPLING_HPP_ */
