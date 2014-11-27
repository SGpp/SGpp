/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#include "NaiveSampleGenerator.hpp"

#include "Random.hpp"

using namespace sg::base;

namespace sg {
  namespace mcm {

      void NaiveSampleGenerator::getSample(sg::base::DataVector& dv) {
        // generate random sample with dimensionality corresponding to the 
        // size of the given datavector (in 0 to 1)
	for(size_t i = 0; i < dv.getSize(); i++)
	{
	  dv[i] = Random::random_double();
	}
      }

    }
  }
