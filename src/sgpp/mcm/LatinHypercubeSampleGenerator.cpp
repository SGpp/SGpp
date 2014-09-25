/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#include "LatinHypercubeSampleGenerator.hpp"

#include <cmath>
#include <list>
#include <algorithm>
#include <vector>
#include "Random.hpp"

using namespace sg::base;

namespace sg {
  namespace mcm {
    
    LatinHypercubeSampleGenerator::LatinHypercubeSampleGenerator(size_t dimensions, size_t numberOfSamples) : SampleGenerator(dimensions) {
      
      // each dimension is divided in n strata to provide n sample points
      // according to latin hypercube sample distribution
      this->numberOfStrata = numberOfSamples;
      this->numberOfSamples = numberOfSamples;
      
      // index number of current sample [1, n]
      this->numberOfCurrentSample = 1;
      
      // equidistant split of [0,1] in n strata -> size of one stratum = 1 / n
      this->sizeOfStrata = 1./ (double) numberOfStrata;
      
      for (size_t i = 0; i < dimensions; i++) {
	currentStrata.push_back(std::vector<size_t>());
	for( size_t j = 0; j < numberOfStrata; j++ ) {
	  currentStrata[i].push_back(j);
	}
      }
      
      shuffleStrataSequence();
    }
    
    void LatinHypercubeSampleGenerator::getSample(sg::base::DataVector& dv) {
      
      // compute random value inside the current stratum selected from the shuffled strata sequence
      for( size_t i = 0; i < dimensions; i++ ) {
	dv[i] = ((double)currentStrata[i][numberOfCurrentSample-1] + Random::random_double()) * sizeOfStrata;
      }
      
      // select next sample from strata sequence. If one sequence is complete shuffle strata to get a new one.
      if(numberOfCurrentSample < numberOfStrata) {
	numberOfCurrentSample++;
      } else {
	numberOfCurrentSample = 0;
	shuffleStrataSequence();
      }
    }
    
    static int random_number_generator(int max) {
        return Random::random() % max;
    }

    void LatinHypercubeSampleGenerator::shuffleStrataSequence() {
      for (size_t i = 0; i < dimensions; i++) {
	std::random_shuffle(currentStrata[i].begin(), currentStrata[i].end(), random_number_generator);
      }
    }    
    
  }
}
