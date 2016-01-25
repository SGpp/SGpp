// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HALTONSAMPLEGENERATOR_HPP
#define HALTONSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/globaldef.hpp>
#include "SampleGenerator.hpp"

namespace SGPP {
  namespace quadrature {

    /**
     *
     */
    class HaltonSampleGenerator : public SampleGenerator {

      public:

        /**
         * Standard constructor
         *
         * @param dimension number of dimensions used for sample generation
         */
        HaltonSampleGenerator(size_t dimension);

        /**
         * Destructor
         */
        ~HaltonSampleGenerator();

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */
        virtual void getSample(SGPP::base::DataVector& sample);

      private:
        int index;
        std::vector<size_t> baseVector;
        std::vector<float_t> iVector;
        std::vector<float_t> fVector;
        std::vector<float_t> resultVector;

    };

  }
}

#endif /* HALTONSAMPLEGENERATOR_HPP */
