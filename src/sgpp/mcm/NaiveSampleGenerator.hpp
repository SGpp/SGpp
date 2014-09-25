/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#ifndef NAIVESAMPLEGENERATOR_HPP
#define NAIVESAMPLEGENERATOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "mcm/SampleGenerator.hpp"

namespace sg {
  namespace mcm {

    /**
     * The class NaiveSampleGenerator implements a simple MonteCarlo sample
     * generator. A sample is generated using the standard random number
     * generator from cmath and transforming the values to double range 0.0 to
     * 1.0.
     */
    class NaiveSampleGenerator : public SampleGenerator {

      public:

        /**
         * Standard constructor 
         *
         * @param dimensions number of dimensions used for sample generation
         */
        NaiveSampleGenerator(size_t dimension): SampleGenerator(dimension){};

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */
        virtual void getSample(sg::base::DataVector& sample);

    };

  }
}

#endif /* NAIVESAMPLEGENERATOR_HPP */
