/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Benjamin Peherstorfer (pehersto@in.tum.de)

#ifndef OPERATIONDENSITYSAMPLING1D_HPP
#define OPERATIONDENSITYSAMPLING1D_HPP

#include "base/grid/Grid.hpp"
#include <cstring>

namespace sg
{
namespace datadriven
{

  /**
   * Sample 1D Probability Density Function
   */
  
  class OperationDensitySampling1D
  {
  public:
    OperationDensitySampling1D() {}
    virtual ~OperationDensitySampling1D() {}

    /**
     * Sample 1D (Density) Functions
     * 
     * @param alpha Coefficient vector for current grid
     * @param mg Referenz of grid pointer 
     * @param malpha Coefficient vector for new grid (mg). Will be resized.
     * @param mdim Marginalize in dimension mdim 
     */
    virtual void doSampling1D(base::DataVector* alpha, size_t num_samples, base::DataVector* samples) = 0;
  };
    
}
}
#endif /* OPERATIONDENSITYSAMPLING1D_HPP */
