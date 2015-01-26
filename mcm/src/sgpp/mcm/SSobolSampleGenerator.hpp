/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#ifndef SSOBOLGENERATOR_HPP
#define SSOBOLGENERATOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "mcm/SampleGenerator.hpp"
#include "mcm/ssobol.h"
#include "Random.hpp"

namespace sg {
  namespace mcm {

    /**
     * A Sample generator based on the scrambled Sobol sequences implemented in SSOBOL.F
     * (ACM Algorithms 823 and 659).
     * Uses a C port of SSOBOL.F by Michael Baudin.
     */
    class SSobolSampleGenerator : public SampleGenerator {

    private:
        Ssobol gen;
        int ok;
    public:

        /**
         * SSOBOL initialisation Status.
         * @return 1 if parameters were ok, else 0.
         */
         int isOk();
          
        /**
         * Standard constructor 
         *
         * @param dimension number of dimensions used for sample generation
         * @param count maximum number of Samples to generate.
         * @param iflag type of scrambling to use, according to SSOBOL:
         *      iflag = 0 : No Scrambling
         *      iflag = 1 : Owen type Scrambling
         *      iflag = 2 : Faure-Tezuka type Scrambling
         *      iflag = 3 : Owen + Faure-Tezuka type Scrambling
         * @param defaultseed Use the (fixed!) seed that is provided with SSOBOL. Default is using a random seed.
         */
        // 30 is recommended in ssobol.h
        SSobolSampleGenerator(size_t dimension, int count, int iflag = 1, bool defaultseed = 0) : SampleGenerator(dimension), gen((int) dimension, count, iflag, 30, &this->ok){
            if (!defaultseed) {
                double seed[24];
                for (int i = 0; i < 24; i++) 
                    seed[i] = Random::random_double();
                // TODO: do not initialize gen twice.
                gen = Ssobol((int) dimension, count, iflag, 30, seed, &this->ok);
            } 
        };

        /**
         * This method generates one sample .
         * Implementation of the abstract Method getSample from SampleGenerator.
         *
         * @param sample DataVector storing the new generated sample vector.
         */
        virtual void getSample(sg::base::DataVector& sample);

    };

  }
}

#endif /* SSOBOLGENERATOR_HPP */
