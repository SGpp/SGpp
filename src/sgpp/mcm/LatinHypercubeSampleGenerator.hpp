/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#ifndef LATINHYPERCUBESAMPLEGENERATOR_HPP
#define LATINHYPERCUBESAMPLEGENERATOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "mcm/SampleGenerator.hpp"

#include <list>
#include <vector>

namespace sg {
  namespace mcm {
    
    /**
     * The class NaiveSampleGenerator implements a simple MonteCarlo sample
     * generator. A sample is generated using the standard random number
     * generator from cmath and transforming the values to double range 0.0 to
     * 1.0.
     */
    class LatinHypercubeSampleGenerator : public SampleGenerator {

      public:

        /**
         * Standard constructor 
         *
         * @param dimensions number of dimensions used for sample generation
         * @param numberOfSamples number of samples to be drawn
         */

        LatinHypercubeSampleGenerator(size_t dimensions, size_t numberOfSamples);

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */

        void getSample(sg::base::DataVector& sample);
        
      private:

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         */

        void shuffleStrataSequence();

        //
        size_t numberOfStrata;

        //       
        size_t numberOfCurrentSample;

        //
        size_t numberOfSamples;
        
        //
        double sizeOfStrata;
  
        //      
        std::vector< std::vector<size_t> > currentStrata;
      
    };
    
  }
}

#endif /* LATINHYPERCUBESAMPLEGENERATOR_HPP */
