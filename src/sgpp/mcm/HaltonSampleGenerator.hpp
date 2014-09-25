/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#ifndef HALTONSAMPLEGENERATOR_HPP
#define HALTONSAMPLEGENERATOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "mcm/SampleGenerator.hpp"

namespace sg {
  namespace mcm {

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
        virtual void getSample(sg::base::DataVector& sample);

      private:
	int index;
        int* base_vektor;
	double* i_vektor;
	double* f_vektor;
	double* result_vektor;

    };

  }
}

#endif /* HALTONSAMPLEGENERATOR_HPP */
