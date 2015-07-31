// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STRATIFIEDSAMPLEGENERATOR_HPP
#define STRATIFIEDSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>
#include "SampleGenerator.hpp"

namespace SGPP {
  namespace quadrature {

    /**
     * The class StratifiedSampleGenerator subdivides every dimension in a given
     * number of strata. For each strata one sample point is generated. In case
     * one sample has already been generated for every strata, the next requested
     * sample will be placed into the first strata.
     */
    class StratifiedSampleGenerator: public SampleGenerator {

      public:

        /**
         * Standard constructor
         *
         * @param dimensions number of dimensions used for sample generation
         * @param strataPerDimension array holding the number of strata used to
         * subdivide the specific dimension
         * @param seed seed for random number generator; if it is equal to -1 the current time is taken as seed
         */

        StratifiedSampleGenerator(std::vector<size_t>& strataPerDimension,
                                  int seed = -1);

        /**
         * Destructor
         */
        virtual ~StratifiedSampleGenerator();

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */

        void getSample(SGPP::base::DataVector& sample);

      private:
        // Array containing the number of strata per dimension
        std::vector<size_t> numberOfStrata;
        // Array containing the current strata number for every dimension
        std::vector<size_t> currentStrata;

        // total number of samples which can be generated for given strata configuration
        size_t numberOfSamples;
        // index number of current sample [0..numberOfSamples-1]
        size_t numberOfCurrentSample;

        // Array containing the size of dimension i strata when dividing [0,1] into numberOfStrata[i]
        std::vector<float_t> sizeOfStrata;

        /**
         * This method computes in which strata the next sample should be generated.
         * Dimension after dimension each stratum is used to generate one sample point. As soon
         * as one dimension is completed the algorithm will start at the beginning of this dimension and
         * counts up the next dimension by 1.
         */
        void getNextStrata();

    };

  }
}

#endif /* STRATIFIEDSAMPLEGENERATOR_HPP */
