/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#ifndef HALTONSAMPLEGENERATOR_HPP
#define HALTONSAMPLEGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/mcm/SampleGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        virtual void getSample(SGPP::base::DataVector& sample);

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
