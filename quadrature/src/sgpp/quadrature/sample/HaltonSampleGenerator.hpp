// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef HALTONSAMPLEGENERATOR_HPP
#define HALTONSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sample/SampleGenerator.hpp>


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
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampelGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */
        virtual void getSample(SGPP::base::DataVector& sample);

      private:
	int index;
        int* base_vektor;
	float_t* i_vektor;
	float_t* f_vektor;
	float_t* result_vektor;

    };

  }
}

#endif /* HALTONSAMPLEGENERATOR_HPP */