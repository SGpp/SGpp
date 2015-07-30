// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LATINHYPERCUBESAMPLEGENERATOR_HPP
#define LATINHYPERCUBESAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <list>
#include <vector>

#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sample/SampleGenerator.hpp>


namespace SGPP {
  namespace quadrature {

    /**
     * The class NaiveSampleGenerator implements a simple MonteCarlo sample
     * generator. A sample is generated using the standard random number
     * generator from cmath and transforming the values to float_t range 0.0 to
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

        LatinHypercubeSampleGenerator(size_t dimensions, size_t numberOfSamples,
                int seed = -1);

        /**
         * Destructor
         */
        virtual ~LatinHypercubeSampleGenerator();

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */

        void getSample(SGPP::base::DataVector& sample);

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
        float_t sizeOfStrata;

        //
        std::vector< std::vector<size_t> > currentStrata;

    };

  }
}

#endif /* LATINHYPERCUBESAMPLEGENERATOR_HPP */
