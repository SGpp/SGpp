// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sample/StratifiedSampleGenerator.hpp>

#include <cmath>

#include <sgpp/quadrature/Random.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace quadrature {
    
    StratifiedSampleGenerator::StratifiedSampleGenerator(size_t dimensions, long long int* strataPerDimension) : SampleGenerator(dimensions)
    {
      
      //TODO: Check correct cast of size_t[dim] and n[i]
      numberOfStrata = ((size_t*)new size_t[dimensions]);
      for(size_t i = 0; i < dimensions; i++) {
	numberOfStrata[i] = ((size_t)strataPerDimension[i]);
      }
      
      // initialise local variables
      currentStrata = new size_t[dimensions];
      sizeOfStrata = new double[dimensions];
      
      // set counter to the first strata for each dimension
      // compute size of strata per dimension
      for(size_t i = 0; i < dimensions; i++) {
	if(numberOfStrata[i] < 1) numberOfStrata[i] = 1;
	currentStrata[i] = 0;
	sizeOfStrata[i] = 1./static_cast<double>(this->numberOfStrata[i]);
      }
      
    }
    
    void StratifiedSampleGenerator::getSample(SGPP::base::DataVector& dv) {
      
      // Check for correct dimension of the parameter vector
      if( dv.getSize() != dimensions ) return;
      
      // Choose a random number inside the stratum selected for this dimension
      for(size_t i = 0; i < dimensions; i++)
      {
	dv[i] = ((double)currentStrata[i] + Random::random_double() ) * sizeOfStrata[i];
      }
      
      // continue to the next stratum used for the next sample
      getNextStrata();
    }
    
    void StratifiedSampleGenerator::getNextStrata() {
      for(size_t i = 0; i < dimensions; i++) {
        // next stratum in this dimension available
	if(currentStrata[i] < numberOfStrata[i] - 1) {
	  currentStrata[i]++;
	  break;
	} else {
          // no more stratum in this dimension, start at the first stratum of this
          // dimension and go on with iterating over the next dimension
	  currentStrata[i] = 0;
	}
      }
    }
    
  }
}